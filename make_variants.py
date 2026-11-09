#!/usr/bin/env python3
"""
make_variants.py — derive model variants of the benchmark XMLs

From each benchmark-{type}-{number}-min.xml (BEAGLE) and -mat.xml (MAT) this writes, beside them:

  mat      benchmark-{type}-{number}-mat.xml
           The BEAGLE tree likelihood replaced by the MAT (mutation-annotated tree) likelihood, which
           keeps an explicit mutation history. The history is resampled along with the SPR moves and
           by a mutationHistoryGibbsOperator at the fixedHeightSPR weight. Derived from -min only.

  gamma    benchmark-{type}-{number}-gamma.xml / -gamma-mat.xml
           A discrete gamma site rate model with 4 categories. The shape (alpha) has an exponential
           prior (mean 0.5) and a scale operator. The MAT version also samples each variable site's rate
           category with a siteRateCategoryGibbsOperator at 1% of the SubtreeLeap weight (the categories
           mix quickly; see BEAST optimization_notes/mat_site_rate_categories/DESIGN.md).

  skygrid  benchmark-{type}-{number}-skygrid.xml / -skygrid-mat.xml
           The constant-size coalescent replaced by a Bayesian skygrid. The constant-size model is kept
           only to simulate the starting tree, so the starting state is that of the source XML.
           Grid intervals: about one per 20 coalescent events (taxa - 1), at least 10 and at most 100.
           Cut-off (grid depth): see CUT_OFFS.

The header block is rewritten to describe the variant. Initial Joint density is filled in by running
BEAST (with --beast) to state 0 at the header's seed. Final Joint density is left out until a full run
of each variant has established it.

usage: make_variants.py [--variants mat,gamma,skygrid] [--beast "java ... -jar beast.jar"] [--datasets small-1,...]
"""
import argparse
import math
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARKS  = _SCRIPT_DIR / "Benchmarks"

# Skygrid cut-off (time before the most recent tip, in the tree's units) for each dataset: the 97.5%
# quantile of the root height pooled over the M5 Max runs of 2026-10-18/19 (BEAGLE -min and MAT -mat,
# 10% burn-in), rounded up to 2 significant figures. Those quantiles were:
#   small-1 1152 (BEAGLE) / 1134 (MAT)      small-2 7.53 / 7.55       tall-1 8.73 / 8.53
#   big-1   0.873 / 0.794                   big-2   2.87 / 2.23
#   wide-1  0.0200 / 0.0200 (undated tips: substitutions per site)
CUT_OFFS = {
    "small-1": 1200.0,
    "small-2": 7.6,
    "tall-1":  8.8,
    "big-1":   0.88,
    "big-2":   2.9,
    "wide-1":  0.020,
    # big-3 and tall-2 have never been run, so these two are estimates rather than measurements.
    # big-3 is SARS-CoV-2 sampled 2020.08-2020.83, whose root sits near the start of the pandemic,
    # about 0.9 years before the most recent tip.
    "big-3":   0.95,
    # tall-2's sampling span is 12.3 years (2000.0-2012.3) and tall-1's root sat about 3.5 years beyond its oldest tip.
    # Replace it with the measured quantile once tall-2 has a chain to take it from.
    "tall-2":  16.0,
}


def grid_intervals(taxa: int) -> int:
    """About one skygrid interval per 20 coalescent events, between 10 and 100."""
    return max(10, min(100, round((taxa - 1) / 20)))


def replace_once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    if n != 1:
        raise ValueError(f"expected one occurrence of {old!r}, found {n}")
    return text.replace(old, new)


def header_field(text: str, field: str) -> str:
    return re.search(rf"^{field}:\s*(.*?)\s*$", text, re.MULTILINE).group(1)


def set_header_field(text: str, field: str, value: str) -> str:
    return re.sub(rf"^{field}:.*$", lambda m: f"{field}: {value}", text, count=1, flags=re.MULTILINE)


def add_header_operators(text: str, lines: list) -> str:
    return re.sub(r"^(MCMC length:)", lambda m: "".join(f"  - {l}\n" for l in lines) + m.group(1),
                  text, count=1, flags=re.MULTILINE)


def rename_outputs(text: str, name: str) -> str:
    text = re.sub(r'fileName="benchmark-[\w-]+\.log"', f'fileName="{name}.log"', text)
    text = re.sub(r'fileName="benchmark-[\w-]+\.trees"', f'fileName="{name}.trees"', text)
    return re.sub(r'operatorAnalysis="benchmark-[\w-]+\.ops"', f'operatorAnalysis="{name}.ops"', text)


def add_mat(xml: str) -> str:
    """Swap the BEAGLE tree likelihood for the MAT one, and sample the mutation history alongside."""
    xml = replace_once(xml, '<patterns id="patterns" from="1" strip="false">',
                       '<patterns id="patterns" from="1" strip="false" constantOnly="true">')

    likelihood = re.search(r'\t<!-- Likelihood for tree given sequence data[^\n]*\n'
                           r'\t<treeDataLikelihood id="treeLikelihood".*?</treeDataLikelihood>\n',
                           xml, re.DOTALL)
    if not likelihood:
        raise ValueError("no treeDataLikelihood block to replace")
    xml = xml.replace(likelihood.group(0), """\t<!-- Mutation-annotated tree (MAT) likelihood for tree given sequence data   -->
	<matDataLikelihood id="treeLikelihood">
		<patterns idref="patterns"/>
		<treeModel idref="treeModel"/>
		<HKYModel idref="jc"/>
		<siteModel idref="siteModel"/>
		<strictClockBranchRates idref="branchRates"/>
	</matDataLikelihood>

	<!-- Resamples the ancestral states (mutation history) along with SPR moves   -->
	<mutationHistoryResampler id="matResampler">
		<matDataLikelihood idref="treeLikelihood"/>
	</mutationHistoryResampler>

	<!-- Number of mutations and ambiguities in the MAT                          -->
	<matStatistic id="mat" name="mat">
		<matDataLikelihood idref="treeLikelihood"/>
	</matStatistic>
""")

    # The Gibbs move over the history runs at the same weight as the SPR it accompanies
    weight = re.search(r'<fixedHeightSubtreePruneRegraft weight="([^"]+)"', xml).group(1)
    xml = replace_once(xml, "\t</operators>", f"""\t\t<mutationHistoryGibbsOperator weight="{weight}" size="3">
			<mutationHistoryResampler idref="matResampler"/>
		</mutationHistoryGibbsOperator>
	</operators>""")

    xml = xml.replace('<treeDataLikelihood idref="treeLikelihood"/>',
                      '<matDataLikelihood idref="treeLikelihood"/>')
    xml = re.sub(r'(<log id="fileLog"[^>]*>\n)', r'\1\t\t\t<matStatistic idref="mat"/>\n', xml, count=1)

    xml = set_header_field(xml, "Model",
                           "MAT (mutation-annotated tree) likelihood, " + header_field(xml, "Model"))
    xml = re.sub(r"^(- (?:SubtreeLeap|fixedHeightSPR) on tree with weight [\d.]+)$",
                 r"\1 (with mutation history resampling)", xml, flags=re.MULTILINE)
    return re.sub(r"^(- fixedHeightSPR on tree[^\n]*)$",
                  lambda m: m.group(1) + f"\n- mutationHistoryGibbs on the mutational history with "
                                         f"weight {weight} (a node and up to 2 ancestors)",
                  xml, count=1, flags=re.MULTILINE)


def add_gamma(xml: str, is_mat: bool) -> str:
    xml = replace_once(xml, "\t\t</substitutionModel>\n\t</siteModel>", """\t\t</substitutionModel>
		<gammaShape gammaCategories="4">
			<parameter id="alpha" value="0.5" lower="0.0"/>
		</gammaShape>
	</siteModel>""")
    xml = replace_once(xml, '\t\t\t\t<coalescentLikelihood idref="coalescent"/>', """\t\t\t\t<exponentialPrior mean="0.5" offset="0.0">
					<parameter idref="alpha"/>
				</exponentialPrior>
				<coalescentLikelihood idref="coalescent"/>""")
    operators = """\t\t<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="alpha"/>
		</scaleOperator>
"""
    header_ops = ["Scale on alpha with weight 3.0"]
    if is_mat:
        leap = float(re.search(r'<subtreeLeap [^>]*weight="([^"]+)"', xml).group(1))
        weight = round(0.01 * leap, 3)
        operators += f"""\t\t<siteRateCategoryGibbsOperator weight="{weight:g}">
			<matDataLikelihood idref="treeLikelihood"/>
		</siteRateCategoryGibbsOperator>
"""
        header_ops.append(f"siteRateCategoryGibbs on the variable sites' rate categories with weight {weight:g}")
    xml = replace_once(xml, "\t</operators>", operators + "\t</operators>")
    xml = replace_once(xml, '\t\t\t<parameter idref="clock.rate"/>\n\t\t\t<rateStatistic idref="meanRate"/>',
                       '\t\t\t<parameter idref="clock.rate"/>\n\t\t\t<parameter idref="alpha"/>\n\t\t\t<rateStatistic idref="meanRate"/>')

    model = header_field(xml, "Model")
    xml = set_header_field(xml, "Model", replace_once(model, ", strict clock", " + discrete gamma (4 categories), strict clock"))
    xml = set_header_field(xml, "Priors", header_field(xml, "Priors") + ", alpha ~ Exponential(mean 0.5)")
    return add_header_operators(xml, header_ops)


def add_skygrid(xml: str, dataset: str, taxa: int) -> str:
    intervals, cut_off = grid_intervals(taxa), CUT_OFFS[dataset]

    xml = replace_once(xml, """	<!-- A prior assumption that the population size has remained constant       -->
	<!-- throughout the time spanned by the genealogy.                           -->
	<constantSize id="constant\"""", """	<!-- A constant-size coalescent, used only to simulate the starting tree     -->
	<constantSize id="constant\"""")

    coalescent = re.search(r'\t<!-- Generate a coalescent likelihood[^\n]*\n\t<coalescentLikelihood id="coalescent">.*?</coalescentLikelihood>\n',
                           xml, re.DOTALL)
    if not coalescent:
        coalescent = re.search(r'\t<coalescentLikelihood id="coalescent">.*?</coalescentLikelihood>\n', xml, re.DOTALL)
    xml = xml.replace(coalescent.group(0), f"""	<!-- Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process       -->
	<!-- {intervals} intervals back to a cut-off of {cut_off:g} before the most recent tip -->
	<gmrfSkyGridLikelihood id="skygrid">
		<populationSizes>
			<!-- skygrid.logPopSize is in log units unlike other popSize                 -->
			<parameter id="skygrid.logPopSize" dimension="{intervals}" value="1.0"/>
		</populationSizes>
		<precisionParameter>
			<parameter id="skygrid.precision" value="0.1" lower="0.0"/>
		</precisionParameter>
		<numGridPoints>
			<parameter id="skygrid.numGridPoints" value="{intervals - 1:.1f}"/>
		</numGridPoints>
		<cutOff>
			<parameter id="skygrid.cutOff" value="{cut_off:g}"/>
		</cutOff>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</gmrfSkyGridLikelihood>
""")
    xml = replace_once(xml, """\t\t<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="constant.popSize"/>
		</scaleOperator>
""", """\t\t<gmrfGridBlockUpdateOperator scaleFactor="1.0" weight="2">
			<gmrfSkyrideLikelihood idref="skygrid"/>
		</gmrfGridBlockUpdateOperator>
		<scaleOperator scaleFactor="0.75" weight="1">
			<parameter idref="skygrid.precision"/>
		</scaleOperator>
""")
    xml = replace_once(xml, """\t\t\t\t<gammaPrior shape="0.001" scale="1000.0" offset="0.0">
					<parameter idref="constant.popSize"/>
				</gammaPrior>
				<coalescentLikelihood idref="coalescent"/>""", """\t\t\t\t<gammaPrior shape="0.001" scale="1000.0" offset="0.0">
					<parameter idref="skygrid.precision"/>
				</gammaPrior>
				<gmrfSkyGridLikelihood idref="skygrid"/>""")
    xml = replace_once(xml, '\t\t\t<parameter idref="constant.popSize"/>\n', """\t\t\t<parameter idref="skygrid.precision"/>
			<parameter idref="skygrid.logPopSize"/>
			<parameter idref="skygrid.cutOff"/>
""")
    xml = replace_once(xml, '\t\t\t<coalescentLikelihood idref="coalescent"/>\n', '\t\t\t<gmrfSkyGridLikelihood idref="skygrid"/>\n')
    assert "constant.popSize" not in xml.split('</constantSize>', 1)[1]

    model = header_field(xml, "Model")
    xml = set_header_field(xml, "Model", replace_once(
        model, "coalescent constant population tree prior",
        f"Bayesian skygrid coalescent tree prior ({intervals} intervals, cut-off {cut_off:g})"))
    xml = set_header_field(xml, "Priors", replace_once(
        header_field(xml, "Priors"), "popSize ~ Gamma(0.001, 1000)", "skygrid precision ~ Gamma(0.001, 1000)"))
    xml = re.sub(r"^\s*- Scale on popSize with weight [\d.]+\n", "", xml, count=1, flags=re.MULTILINE)
    xml = re.sub(r"upDown on rate and popSize", "upDown on rate and node heights", xml, count=1)
    return add_header_operators(xml, ["gmrfGridBlockUpdate on skygrid.logPopSize with weight 2.0",
                                      "Scale on skygrid.precision with weight 1.0"])


def initial_joint(beast: str, xml: str, seed: int) -> float:
    """Run BEAST to state 0 and return the Joint density it reports there."""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "initial.xml"
        path.write_text(re.sub(r'chainLength="\d+"', 'chainLength="1"', xml, count=1), encoding="utf-8")
        out = subprocess.run(shlex.split(beast) + ["-seed", str(seed), "-overwrite", "-citations_off", str(path)],
                             cwd=tmp, capture_output=True, text=True)
        lines = out.stdout.splitlines()
        try:
            columns = next(l for l in lines if l.startswith("state\t")).split("\t")
            row = next(l for l in lines if l.startswith("0\t")).split("\t")
        except StopIteration:
            raise RuntimeError(f"BEAST did not log state 0:\n{out.stdout[-2000:]}\n{out.stderr[-2000:]}")
        return float(row[[c.strip() for c in columns].index("Joint")])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variants", default="gamma,skygrid")
    parser.add_argument("--datasets", default=",".join(CUT_OFFS))
    parser.add_argument("--beast", default=None,
                        help="BEAST command used to fill in Initial Joint density (omit to leave it out)")
    args = parser.parse_args()

    for dataset in args.datasets.split(","):
        btype, number = dataset.split("-")
        folder = BENCHMARKS / f"Benchmark-{btype}"
        for engine in ("min", "mat"):
            path = folder / f"benchmark-{dataset}-{engine}.xml"
            if not path.exists():
                # A new dataset starts as -min alone; make its -mat first, then the rest
                print(f"no {path.relative_to(_SCRIPT_DIR)}, skipping its variants")
                continue
            source = path.read_text(encoding="utf-8")
            taxa = int(header_field(source, "Taxa"))
            seed = int(header_field(source, "Seed"))
            for variant in args.variants.split(","):
                is_mat = engine == "mat"
                if variant == "mat":
                    if is_mat:          # -mat is what this variant produces, not what it starts from
                        continue
                    name = f"benchmark-{dataset}-mat"
                    xml = add_mat(source)
                else:
                    name = f"benchmark-{dataset}-{variant}" + ("-mat" if is_mat else "")
                    xml = add_gamma(source, is_mat) if variant == "gamma" else add_skygrid(source, dataset, taxa)
                xml = rename_outputs(xml, name)
                xml = set_header_field(xml, "Benchmark", name)
                xml = re.sub(r"^Final Joint density:.*\n", "", xml, count=1, flags=re.MULTILINE)
                if args.beast:
                    xml = set_header_field(xml, "Initial Joint density", f"{initial_joint(args.beast, xml, seed):.4f}")
                else:
                    xml = re.sub(r"^Initial Joint density:.*\n", "", xml, count=1, flags=re.MULTILINE)
                out = folder / f"{name}.xml"
                out.write_text(xml, encoding="utf-8")
                print(f"wrote {out.relative_to(_SCRIPT_DIR)}  {header_field(xml, 'Initial Joint density') if args.beast else ''}")


if __name__ == "__main__":
    main()
