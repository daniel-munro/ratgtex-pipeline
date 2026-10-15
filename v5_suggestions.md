**1. Make the release a configuration, rather than a sequence of manual edits.**

Currently, preparing a release involves changing `config.yaml`, copying Pantry into tissue directories, editing those copies’ configs, running setup scripts, assembling merged tissues, and switching analyses on and off. The top-level DAG does not cover that whole process.

I would define v5 in one versioned configuration, with explicit stage targets such as:

- Alignment and identity QC.
- Phenotyping and tissue assembly.
- QTL mapping and TWAS weights.
- Release exports and validation.

Separate stage workflows are fine; they need explicit inputs and outputs. You should be able to request a stage without editing `rule all` or commenting analyses in individual tissue configs.

**2. Separate sequencing datasets, biological tissues, and analysis cohorts.**

These concepts currently share the “tissue” namespace. NAcc1–4 are source datasets; NAcc is the combined analysis cohort. In [steps/setup.smk](/data/hps/home/dmunro/gdml/ratgtex/steps/setup.smk:38), final analysis tissues are partly inferred from whether their names end in a digit.

I would represent those relationships directly:

- Sequencing datasets describe source studies and library properties.
- Samples connect libraries to animals and tissues.
- Analysis cohorts explicitly list their included samples.

That removes naming assumptions and makes new tissues, multiple studies, and alternative sample subsets much easier to handle.

**3. Make the sample tracking system authoritative.**

The existing original/current maps, relabeling decisions, removal reasons, and `sample_tracking.py` are a useful foundation. However, processing still reads separate release-local `fastq_map.txt` and `rat_ids.txt` files.

I would generate those processing inputs from a versioned sample manifest and decision tables. Use stable library/sample identifiers alongside the corrected animal RFID, so correcting biological identity does not erase the original sequencing identity.

This also makes repeated libraries explicit. The current merged-tissue script combines animal lists with `sort | uniq`; that does not establish how multiple libraries from the same animal should contribute to the merged phenotype.

**4. Use one pinned Pantry checkout, with generated per-cohort configuration.**

The setup scripts copy Pantry source into every tissue directory using `rsync`. Those copies can diverge, and the release directories containing them are ignored by Git.

I would keep a pinned Pantry revision in one location and run it with separate configurations and output directories. Record the Pantry revision and any local modifications in the release provenance.

There is already configuration drift: the Pheast template enables QTL mapping and cross-modality analysis, whereas the Brain working config currently enables TWAS. That may reflect legitimate successive runs, but the final config alone cannot reconstruct them.

**5. Declare modalities and supported analyses in one place.**

Modality information is repeated in the [main Snakefile](/data/hps/home/dmunro/gdml/ratgtex/Snakefile:13), phenotyping config, Pheast config, and the cross-modality assembly command.

For v5, I would define each modality once, including its grouping, producer, and supported analyses. Cohorts could then select explicit modality sets.

Cross-modality analysis deserves its own definition: which modalities participate, how features map to genes/groups, and which sample set and covariates apply. Adding a modality can change the combined analysis, so that membership should be recorded in the release.

**6. Replace the merged-tissue setup script with declared workflow steps.**

[setup_merged_tissues.sh](/data/hps/home/dmunro/gdml/ratgtex/scripts/setup/setup_merged_tissues.sh) hard-codes v4, tissue membership, and Pantry’s internal intermediate paths. It also requires constituent datasets to have reached the right stage beforehand.

The underlying approach—reuse sample-level quantification, then assemble and normalize the combined cohort—makes sense. I would preserve that distinction while making the dependencies explicit.

In particular, distinguish reusable sample-level products from cohort-dependent products such as splice clusters, filtering, normalization, phenotype PCs, and QTL results. Adding samples should trigger the latter steps appropriately.

**7. Treat scientific settings as release parameters.**

Several consequential choices are embedded in scripts or copied Pantry rules:

- Five genotype PCs and twenty phenotype PCs in the inspected Pheast copy.
- Genotype MAF/MAC filters.
- FDR thresholds.
- Autosomes 1–20.
- Single-end fragment-length assumptions.
- The Brain-derived LD reference and TWAS SNP list.

These choices may be appropriate. I would expose and record them before deciding whether to change them. For combined studies, I would specifically review how study/batch effects are handled and whether the same covariate policy suits every modality.

I would also document the decision to use genotype-aware STAR alignment before identity QC, since identity corrections then affect alignment dependencies.

**8. Give release packaging its own workflow.**

The current default targets cover selected QC and association outputs, but they do not define a complete release. TWAS products, cross-modality results, metadata, and QC decisions are not all represented there.

A release target should assemble an explicit product inventory, including:

- Download tables and metadata.
- Sample and feature counts for every configured tissue/modality.
- Reference, allele, coordinate, and statistical definitions.
- Software/configuration provenance and checksums.
- A concise QC and exclusion summary.

Validation should check scientific consistency—sample correspondence, identifier uniqueness, group mappings, variant alleles, and expected table relationships—and fail when assumptions are violated.

A few smaller issues also deserve attention during that work: the README references a root `environment.yml` that is absent; its GPU rule names do not fully match the current rules; sex-concordance QC expects an expression file with no producer in the inspected top-level workflow; reference index paths encode read length but not reference identity; and the trans-QTL retry increases host memory despite its comment describing a CUDA memory problem.

My priority would be **the sample/cohort model, one release configuration, and a pinned Pantry integration**. Those changes would remove most of the repeated setup work and make the new tissues and modalities easier to add. I would retain the existing scientific methods initially, then assess methodological changes separately so their effects remain interpretable.

# Todo

The items below reflect the decisions following review of the suggestions above and supersede those suggestions where they differ.

- [x] Organize v5 setup and RatGTEx processing around the main configuration, with explicit targets for alignment/identity QC, project setup, and processing after Pantry. Prepare self-contained Pantry projects and configs, but run phenotyping and Pheast manually per dataset/tissue, allowing concurrent projects and individual troubleshooting. Make RatGTEx stage selection possible without editing `rule all`.
- [ ] Define clear terminology and explicit relationships for source studies, tissue-specific datasets, tissue analyses, libraries/samples, animals, and genotype panels. Replace tissue-name suffix heuristics with configured tissue membership.
- [ ] Generate `fastq_map.txt` and `rat_ids.txt` automatically from the sample tracking information and relabeling/removal decisions. Preserve original sequencing identity alongside corrected animal identity, and explicitly handle multiple libraries from the same animal.
- [ ] Review tracked samples omitted from v4 for possible inclusion in v5; record inclusion/exclusion decisions and reasons in the tracking information.
- [ ] Keep Pantry integration simple: copy Pantry into self-contained per-dataset or per-tissue project directories during setup, then generate their phenotyping and Pheast configs from the main RatGTEx config. Record the source Pantry revision and allow later tissue-specific analyses without silently overwriting their configuration or results.
- [ ] Define modalities and grouping metadata in the main config and use those definitions throughout setup, workflow targets, and generated Pantry configs. Generate cross-modality phenotypes from all configured modalities and process the result like any other grouped modality, including covariates and mapping.
- [ ] Integrate tissue merging with the dataset/tissue model and workflow dependencies. Replace the v4-specific merge setup with configuration-driven assembly that reuses sample-level quantification and regenerates tissue-dependent products when membership changes.
- [ ] Consolidate scattered scientific parameters in the main config and propagate them into Pantry configs and RatGTEx rules. Preserve the current scientific settings; batch effects and covariate methodology have already been examined and do not need reassessment in this work.
- [ ] Document and support the identity-QC review cycle: align and check identity; manually review matches, searching the full rat genotypes when needed; update sample tracking and regenerate FASTQ maps; delete affected BAMs and regenerate alignment and dependent outputs before continuing.
- [ ] Improve documentation and processing summaries, including sample/feature counts across configured tissues and modalities, QC/exclusion decisions, and software/configuration provenance. Keep portal file assembly and its product inventory in `../ratgtex-server-data/`; pipeline products also support the paper and other projects, so do not impose a separate release-output inventory here.
- [x] Remove the obsolete root `environment.yml` setup instructions from the README and clarify which environment is needed for RatGTEx-side steps versus Pantry's own environment specifications.
- [x] Update GPU rule names in the README to correspond to `qtl.smk` and current Pantry rules (already completed by the user).
- [ ] Update sex-concordance QC to consume Pantry-produced expression and declare the corresponding workflow dependency.
- [ ] Simplify file names and directory organization where useful, updating affected paths and documentation together.
- [ ] Organize `samples/`, e.g. move GEO-related files into geo subdirectory.
- [ ] Remove or update `test/`.

Retain the current read-length-specific reference index paths, since tissues share the reference genome and annotation. Retain the current trans-QTL retry behavior.

Implementation note: the first item is implemented through named RatGTEx targets and `scripts/run_stage.py`, with generated self-contained Pantry projects, explicit `datasets`/`tissues`, configured modalities, and separate QTL/TWAS entry points. Pantry phenotyping and Pheast are launched manually in individual project directories; the runner handles setup and RatGTEx work before/after Pantry. See README for handoffs and dry-run commands. Source sample-file generation and full sample/tissue identity redesign remain separate work. Setup writes initial configs; later changes to project configs, membership, or links are manual.

Terminology decision: a dataset is a tissue-specific source collection (NAcc1); a tissue is the analysis combining one or more datasets (NAcc). Genotype VCFs are identified by `genotype_panel`, independently of RNA-seq dataset names. Use study for the originating project and cohort for a defined group of animals; no study/cohort config entity is introduced yet. This terminology supersedes the original suggestions above.

Configuration layout: the current main config and tissue metadata live in `config/config.yaml` and `config/tissue_info.tsv`. `config/datasets.txt` and `config/tissues.txt` are generated from the main config during the single `setup` stage. Input/output paths remain relative to the repository root. Older release files are archived under `v4/`; external consumers will be updated separately.
