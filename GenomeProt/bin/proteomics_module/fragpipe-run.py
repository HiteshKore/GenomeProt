#!/usr/bin/env python3

import argparse, glob, os, subprocess, sys
import datetime, logging

workflow_template = """# Workflow:
database.db-path=

crystalc.run-crystalc=false
database.decoy-tag=rev_
diann.channel-normalization-strategy=0
diann.cmd-opts=
diann.generate-msstats=true
diann.heavy=
diann.library=
diann.light=
diann.medium=
diann.min-site-prob=0.75
diann.mod-tag=
diann.q-value=0.01
diann.quantification-strategy=3
diann.quantification-strategy-2=2
diann.run-dia-nn=
diann.run-dia-plex=false
diann.run-specific-protein-q-value=false
diann.unrelated-runs=false
diann.use-predicted-spectra=false
diatracer.corr-threshold=0.3
diatracer.delta-apex-im=0.01
diatracer.delta-apex-rt=3
diatracer.mass-defect-filter=true
diatracer.mass-defect-offset=0.1
diatracer.rf-max=500
diatracer.run-diatracer=false
diatracer.write-intermediate-files=false
diaumpire.AdjustFragIntensity=true
diaumpire.BoostComplementaryIon=false
diaumpire.CorrThreshold=0
diaumpire.DeltaApex=0.2
diaumpire.ExportPrecursorPeak=false
diaumpire.Q1=true
diaumpire.Q2=true
diaumpire.Q3=true
diaumpire.RFmax=500
diaumpire.RPmax=25
diaumpire.RTOverlap=0.3
diaumpire.SE.EstimateBG=false
diaumpire.SE.IsoPattern=0.3
diaumpire.SE.MS1PPM=10
diaumpire.SE.MS2PPM=20
diaumpire.SE.MS2SN=1.1
diaumpire.SE.MassDefectFilter=true
diaumpire.SE.MassDefectOffset=0.1
diaumpire.SE.NoMissedScan=1
diaumpire.SE.SN=1.1
diaumpire.run-diaumpire=false
fpop.coadaptr.fpop.fpop_masses=
fpop.coadaptr.fpop.run-fpop-coadaptr=false
fpop.fragpipe.fpop.fpop-tmt=false
fpop.fragpipe.fpop.label_control=
fpop.fragpipe.fpop.label_fpop=
fpop.fragpipe.fpop.region_size=1
fpop.fragpipe.fpop.run-fpop=false
fpop.fragpipe.fpop.subtract-control=false
freequant.mz-tol=10
freequant.rt-tol=0.4
freequant.run-freequant=false
ionquant.excludemods=
ionquant.formula=
ionquant.heavy=
ionquant.imtol=0.05
ionquant.intensitymode=0
ionquant.ionfdr=0.01
ionquant.light=
ionquant.locprob=0.75
ionquant.maxlfq=1
ionquant.mbr=0
ionquant.mbrimtol=0.05
ionquant.mbrmincorr=0
ionquant.mbrrttol=1
ionquant.mbrtoprun=10
ionquant.medium=
ionquant.minfreq=0
ionquant.minions=1
ionquant.minisotopes=1
ionquant.minscans=3
ionquant.mztol=10
ionquant.normalization=1
ionquant.peptidefdr=1
ionquant.proteinfdr=1
ionquant.requantify=1
ionquant.rttol=0.4
ionquant.run-ionquant=true
ionquant.tp=0
ionquant.uniqueness=0
ionquant.use-labeling=false
ionquant.use-lfq=true
ionquant.writeindex=0
metaproteomics.cmd-line-opts=
metaproteomics.delta-hyperscore=0.0
metaproteomics.host-name=Homo sapiens
metaproteomics.iterations=3
metaproteomics.min-pept-cnt-per-prot=1
metaproteomics.min-uniq-pept-cnt=3
metaproteomics.min-uniq-pept-cnt-per-prot=1
metaproteomics.qvalue=0.01
metaproteomics.run-metaproteomics=false
msbooster.find-best-im-model=false
msbooster.find-best-rt-model=false
msbooster.find-best-spectra-model=false
msbooster.fragmentation-type=0
msbooster.im-model=DIA-NN
msbooster.koina-url=
msbooster.predict-im=true
msbooster.predict-rt=true
msbooster.predict-spectra=true
msbooster.rt-model=DIA-NN
msbooster.run-msbooster=true
msbooster.spectra-model=DIA-NN
msbooster.spectral-library-path=
msfragger.Y_type_masses=
msfragger.activation_types=all
msfragger.allowed_missed_cleavage_1=1
msfragger.allowed_missed_cleavage_2=1
msfragger.analyzer_types=all
msfragger.calibrate_mass=2
msfragger.check_spectral_files=false
msfragger.clip_nTerm_M=true
msfragger.deisotope=1
msfragger.delta_mass_exclude_ranges=(-1.5,3.5)
msfragger.deneutralloss=1
msfragger.diagnostic_fragments=
msfragger.diagnostic_intensity_filter=0
msfragger.digest_max_length=50
msfragger.digest_min_length=7
msfragger.fragment_ion_series=b,y
msfragger.fragment_mass_tolerance=20
msfragger.fragment_mass_units=1
msfragger.group_variable=0
msfragger.intensity_transform=0
msfragger.ion_series_definitions=
msfragger.isotope_error=0/1/2
msfragger.labile_search_mode=off
msfragger.localize_delta_mass=false
msfragger.mass_diff_to_variable_mod=0
msfragger.mass_offsets=0
msfragger.mass_offsets_detailed=
msfragger.max_fragment_charge=2
msfragger.max_variable_mods_combinations=5000
msfragger.max_variable_mods_per_peptide=3
msfragger.min_fragments_modelling=2
msfragger.min_matched_fragments=4
msfragger.min_sequence_matches=2
msfragger.minimum_peaks=15
msfragger.minimum_ratio=0.00
msfragger.misc.fragger.clear-mz-hi=0
msfragger.misc.fragger.clear-mz-lo=0
msfragger.misc.fragger.digest-mass-hi=5000
msfragger.misc.fragger.digest-mass-lo=500
msfragger.misc.fragger.enzyme-dropdown-1=
msfragger.misc.fragger.enzyme-dropdown-2=null
msfragger.misc.fragger.precursor-charge-hi=4
msfragger.misc.fragger.precursor-charge-lo=2
msfragger.misc.fragger.remove-precursor-range-hi=1.5
msfragger.misc.fragger.remove-precursor-range-lo=-1.5
msfragger.misc.slice-db=1
msfragger.num_enzyme_termini=2
msfragger.output_format=pepXML_pin
msfragger.output_max_expect=50
msfragger.output_report_topN=1
msfragger.output_report_topN_dda_plus=5
msfragger.output_report_topN_dia1=5
msfragger.override_charge=false
msfragger.precursor_mass_lower=-20
msfragger.precursor_mass_mode=selected
msfragger.precursor_mass_units=1
msfragger.precursor_mass_upper=20
msfragger.precursor_true_tolerance=20
msfragger.precursor_true_units=1
msfragger.remainder_fragment_masses=
msfragger.remove_precursor_peak=1
msfragger.report_alternative_proteins=true
msfragger.require_precursor=true
msfragger.restrict_deltamass_to=all
msfragger.reuse_dia_fragment_peaks=false
msfragger.run-msfragger=true
msfragger.search_enzyme_cut_1=
msfragger.search_enzyme_cut_2=
msfragger.search_enzyme_name_1=
msfragger.search_enzyme_name_2=null
msfragger.search_enzyme_nocut_1=
msfragger.search_enzyme_nocut_2=
msfragger.search_enzyme_sense_1=
msfragger.search_enzyme_sense_2=C
msfragger.table.fix-mods=0.0,C-Term Peptide,true,-1; 0.0,N-Term Peptide,true,-1; 0.0,C-Term Protein,true,-1; 0.0,N-Term Protein,true,-1; 0.0,G (glycine),true,-1; 0.0,A (alanine),true,-1; 0.0,S (serine),true,-1; 0.0,P (proline),true,-1; 0.0,V (valine),true,-1; 0.0,T (threonine),true,-1; 57.02146,C (cysteine),true,-1; 0.0,L (leucine),true,-1; 0.0,I (isoleucine),true,-1; 0.0,N (asparagine),true,-1; 0.0,D (aspartic acid),true,-1; 0.0,Q (glutamine),true,-1; 0.0,K (lysine),true,-1; 0.0,E (glutamic acid),true,-1; 0.0,M (methionine),true,-1; 0.0,H (histidine),true,-1; 0.0,F (phenylalanine),true,-1; 0.0,R (arginine),true,-1; 0.0,Y (tyrosine),true,-1; 0.0,W (tryptophan),true,-1; 0.0,B ,true,-1; 0.0,J,true,-1; 0.0,O,true,-1; 0.0,U,true,-1; 0.0,X,true,-1; 0.0,Z,true,-1
msfragger.table.var-mods=15.9949,M,true,1; 42.0106,[^,true,1; 79.96633,STY,false,3; -17.0265,nQnC,false,1; -18.0106,nE,false,1; 0.0,site_06,false,1; 0.0,site_07,false,1; 0.0,site_08,false,1; 0.0,site_09,false,1; 0.0,site_10,false,1; 0.0,site_11,false,1; 0.0,site_12,false,1; 0.0,site_13,false,1; 0.0,site_14,false,1; 0.0,site_15,false,1; 0.0,site_16,false,1
msfragger.track_zero_topN=0
msfragger.use_all_mods_in_first_search=false
msfragger.use_detailed_offsets=false
msfragger.use_topN_peaks=1000
msfragger.write_calibrated_mzml=false
msfragger.zero_bin_accept_expect=0
msfragger.zero_bin_mult_expect=1
opair.activation1=HCD
opair.activation2=ETD
opair.allowed_sites=
opair.filterOxonium=true
opair.glyco_db=
opair.max_glycans=4
opair.max_isotope_error=2
opair.min_isotope_error=0
opair.ms1_tol=20
opair.ms2_tol=20
opair.oxonium_filtering_file=
opair.oxonium_minimum_intensity=0.05
opair.reverse_scan_order=false
opair.run-opair=false
opair.single_scan_type=false
peptide-prophet.cmd-opts=--decoyprobs --ppm --accmass --nonparam --expectscore
peptide-prophet.combine-pepxml=false
peptide-prophet.run-peptide-prophet=false
percolator.cmd-opts=--no-terminate --post-processing-tdc --subset-max-train 500000
percolator.keep-tsv-files=false
percolator.min-prob=0.7
percolator.run-percolator=true
phi-report.dont-use-prot-proph-file=false
phi-report.filter=--picked --prot 0.01 --minPepLen 8
phi-report.pep-level-summary=false
phi-report.print-decoys=false
phi-report.prot-level-summary=false
phi-report.remove-contaminants=false
phi-report.run-report=true
protein-prophet.cmd-opts=--maxppmdiff 2000000 --minprob 0.5
protein-prophet.run-protein-prophet=true
ptmprophet.cmdline=
ptmprophet.override-defaults=false
ptmprophet.run-ptmprophet=false
ptmshepherd.adv_params=false
ptmshepherd.annotate_assigned_mods=false
ptmshepherd.annotation-common=false
ptmshepherd.annotation-custom=false
ptmshepherd.annotation-glyco=false
ptmshepherd.annotation-unimod=true
ptmshepherd.annotation_file=
ptmshepherd.annotation_tol=0.01
ptmshepherd.cap_y_ions=
ptmshepherd.decoy_type=1
ptmshepherd.diag_ions=
ptmshepherd.diagmine_diagMinFoldChange=3.0
ptmshepherd.diagmine_diagMinSpecDiff=25
ptmshepherd.diagmine_fragMinFoldChange=3.0
ptmshepherd.diagmine_fragMinPropensity=12.5
ptmshepherd.diagmine_fragMinSpecDiff=25
ptmshepherd.diagmine_minIonsPerSpec=2
ptmshepherd.diagmine_minPeps=25
ptmshepherd.diagmine_pepMinFoldChange=3.0
ptmshepherd.diagmine_pepMinSpecDiff=25
ptmshepherd.glyco_fdr=1.00
ptmshepherd.glyco_isotope_max=3
ptmshepherd.glyco_isotope_min=-1
ptmshepherd.glyco_ppm_tol=50
ptmshepherd.glycodatabase=
ptmshepherd.histo_smoothbins=2
ptmshepherd.iontype_a=false
ptmshepherd.iontype_b=true
ptmshepherd.iontype_c=true
ptmshepherd.iontype_x=false
ptmshepherd.iontype_y=true
ptmshepherd.iontype_z=true
ptmshepherd.localization_allowed_res=
ptmshepherd.n_glyco=true
ptmshepherd.normalization-psms=true
ptmshepherd.normalization-scans=false
ptmshepherd.output_extended=false
ptmshepherd.peakpicking_mass_units=0
ptmshepherd.peakpicking_minPsm=10
ptmshepherd.peakpicking_promRatio=0.3
ptmshepherd.peakpicking_width=0.002
ptmshepherd.precursor_mass_units=0
ptmshepherd.precursor_tol=0.01
ptmshepherd.print_decoys=false
ptmshepherd.print_full_glyco_params=false
ptmshepherd.prob_mass=0.5
ptmshepherd.remainder_masses=
ptmshepherd.remove_glycan_delta_mass=true
ptmshepherd.run-shepherd=false
ptmshepherd.run_diagextract_mode=false
ptmshepherd.run_diagmine_mode=false
ptmshepherd.run_glyco_mode=false
ptmshepherd.spectra_condPeaks=150
ptmshepherd.spectra_condRatio=0.0001
ptmshepherd.spectra_maxPrecursorCharge=4
ptmshepherd.spectra_maxfragcharge=2
ptmshepherd.spectra_ppmtol=20
ptmshepherd.use_msfragger_localization=false
ptmshepherd.varmod_masses=
quantitation.run-label-free-quant=false
run-psm-validation=true
run-validation-tab=true
saintexpress.cmd-opts=
saintexpress.max-replicates=10
saintexpress.run-saint-express=false
saintexpress.virtual-controls=100
skyline.run-skyline=false
skyline.skyline=false
skyline.skyline-custom=false
skyline.skyline-custom-path=
skyline.skyline-daily=true
skyline.skyline-fragment-tolerance=10
skyline.skyline-mods-mode=Default
skyline.skyline-precursor-tolerance=10
skyline.use-ssl=false
speclibgen.convert-pepxml=true
speclibgen.convert-psm=false
speclibgen.easypqp.extras.max_delta_ppm=15
speclibgen.easypqp.extras.max_delta_unimod=0.02
speclibgen.easypqp.extras.max_glycan_qval=1
speclibgen.easypqp.extras.rt_lowess_fraction=0
speclibgen.easypqp.fragment.a=false
speclibgen.easypqp.fragment.b=true
speclibgen.easypqp.fragment.c=false
speclibgen.easypqp.fragment.x=false
speclibgen.easypqp.fragment.y=true
speclibgen.easypqp.fragment.z=false
speclibgen.easypqp.im-cal=Automatic selection of a run as reference IM
speclibgen.easypqp.labile_mode=Regular (not glyco)
speclibgen.easypqp.neutral_loss=false
speclibgen.easypqp.rt-cal=noiRT
speclibgen.easypqp.select-file.text=
speclibgen.easypqp.select-im-file.text=
speclibgen.keep-intermediate-files=false
speclibgen.run-speclibgen=true
tab-run.delete_temp_files=false
tab-run.export_matched_fragments=false
tab-run.sub_mzml_prob_threshold=0.5
tab-run.write_sub_mzml=false
tmtintegrator.add_Ref=-1
tmtintegrator.aggregation_method=0
tmtintegrator.allow_overlabel=true
tmtintegrator.allow_unlabeled=true
tmtintegrator.best_psm=true
tmtintegrator.channel_num=TMT-6
tmtintegrator.extraction_tool=IonQuant
tmtintegrator.glyco_qval=-1
tmtintegrator.groupby=-1
tmtintegrator.log2transformed=true
tmtintegrator.max_pep_prob_thres=0
tmtintegrator.min_ntt=0
tmtintegrator.min_pep_prob=0.9
tmtintegrator.min_percent=0.05
tmtintegrator.min_purity=0.5
tmtintegrator.min_resolution=0
tmtintegrator.min_site_prob=-1
tmtintegrator.min_snr=0
tmtintegrator.mod_tag=none
tmtintegrator.ms1_int=true
tmtintegrator.outlier_removal=true
tmtintegrator.philosopher-msstats=false
tmtintegrator.print_RefInt=false
tmtintegrator.prot_exclude=none
tmtintegrator.prot_norm=0
tmtintegrator.psm_norm=false
tmtintegrator.quant_level=2
tmtintegrator.ref_d_tag=Pool
tmtintegrator.ref_tag=Bridge
tmtintegrator.run-tmtintegrator=false
tmtintegrator.tolerance=20
tmtintegrator.unique_gene=0
tmtintegrator.unique_pep=false
tmtintegrator.use_glycan_composition=false
workflow.description=<p>Custom FragPipe workflow generated from the GenomeProt proteomics module.</p>
workflow.input.data-type.im-ms=false
workflow.input.data-type.regular-ms=true
workflow.misc.save-sdrf=true
workflow.misc.sdrf-type=Default
workflow.saved-with-ver=23.1
""".replace('\r', '')

protease_rules = {
    "stricttrypsin": ("KR", '', 'C'),
    "trypsin": ("KR", 'P', 'C'),
    "trypsin_gluc": ("DEKR", 'P', 'C'),
    "gluc": ("DE", 'P', 'C'),
    "lysc": ('K', 'P', 'C'),
    "lysn": ('K', '', 'N'),
    "argc": ('R', 'P', 'C'),
    "aspn": ('D', '', 'N')
}

manifest_path = "data.fp-manifest"
workflow_path = ''

accepted_mass_spec_data_types = ("DDA", "DDA+", "DIA", "DIA-Quant", "DIA-Lib", "GPF-DIA")
accepted_mass_spec_data_types_lower = tuple(accepted_mass_spec_data_type.lower() for accepted_mass_spec_data_type in accepted_mass_spec_data_types)

# Run a command and log all stdout and stderr messages into a log file
def run_command(args):
    with subprocess.Popen(args, text = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT) as proc:
        for raw_line in proc.stdout:
            line = raw_line.strip()
            if line:
                logging.info(line)

# Turn a relative path into a real path whenever possible
def get_real_path(path):
    real_path = path
    try:
        real_path = os.path.realpath(path)
    except:
        pass
    return real_path

# Add decoys (and contaminants) to a proteome database for the peptide search
def prepare_proteome_database(db_path, philosopher_path, is_add_contaminants):
    try:
        logging.info("Initializing the workspace...")
        run_command([philosopher_path, "workspace", "--init"])
        logging.info(f"Preparing the custom proteome database by adding decoys{' and contaminants' if is_add_contaminants else ''}...")
        run_command([philosopher_path, "database", "--custom", db_path] + ["--contam"] * is_add_contaminants)
        custom_db_filenames = glob.glob("*.fas")
        if len(custom_db_filenames) == 0:
            logging.error("Error: The custom proteome database could not be prepared.")
            return False
        custom_db_filename = custom_db_filenames[0]
        os.rename(custom_db_filename, "custom_database.fasta")
        logging.info("Cleaning up the workspace...")
        run_command([philosopher_path, "workspace", "--clean"])
        logging.info("Custom proteome database prepared at 'custom_database.fasta'.")
        return True
    except Exception as error:
        logging.error("Error: An exception occurred when preparing the custom proteome database.")
        logging.error(error)
        return False

# Generate the manifest file for the peptide search
def prepare_manifest(mass_spec_filenames, mass_spec_data_types):
    try:
        logging.info("Generating the manifest file...")
        data = '\n'.join([f"{mass_spec_filename}\t\t\t{mass_spec_data_type}" for (mass_spec_filename, mass_spec_data_type) in zip(mass_spec_filenames, mass_spec_data_types)]) + '\n'
        with open(manifest_path, 'w') as f:
            _ = f.write(data)
        logging.info(f"Manifest file generated at '{manifest_path}'.")
        return True
    except Exception as error:
        logging.error("Error: An exception occurred when generating the manifest file.")
        logging.error(error)
        return False

# Generate the workflow file for the peptide search
def prepare_workflow(protease1, protease2, is_perform_quantification):
    global workflow_path
    try:
        logging.info("Generating the workflow file...")
        custom_db_path = os.path.join(os.getcwd(), "custom_database.fasta")
        if not os.path.exists(custom_db_path):
            logging.error(f"Error: The custom proteome database '{custom_db_path}' does not exist.")
            return False
        if not os.path.isfile(custom_db_path):
            logging.error(f"Error: The custom proteome database '{custom_db_path}' exists but it is not a file.")
            return False
        if os.path.getsize(custom_db_path) == 0:
            logging.error(f"Error: The custom proteome database '{custom_db_path}' is empty.")
            return False
        # Specify the proteome database
        workflow_data = workflow_template.replace("database.db-path=", "database.db-path=" + custom_db_path)
        # Specify the proteases
        (cuts, nocuts, sense) = protease_rules[protease1]
        workflow_data = workflow_data.replace("Workflow:", f"Workflow: {protease1}") \
                                     .replace("msfragger.misc.fragger.enzyme-dropdown-1=", "msfragger.misc.fragger.enzyme-dropdown-1=" + protease1) \
                                     .replace("msfragger.search_enzyme_name_1=", "msfragger.search_enzyme_name_1=" + protease1) \
                                     .replace("msfragger.search_enzyme_cut_1=", "msfragger.search_enzyme_cut_1=" + cuts) \
                                     .replace("msfragger.search_enzyme_nocut_1=", "msfragger.search_enzyme_nocut_1=" + nocuts) \
                                     .replace("msfragger.search_enzyme_sense_1=", "msfragger.search_enzyme_sense_1=" + sense)
        if protease2:
            (cuts, nocuts, sense) = protease_rules[protease2]
            workflow_data = workflow_data.replace(f"Workflow: {protease1}", f"Workflow: {protease1}_{protease2}") \
                                         .replace("msfragger.misc.fragger.enzyme-dropdown-2=null", "msfragger.misc.fragger.enzyme-dropdown-2=" + protease2) \
                                         .replace("msfragger.search_enzyme_name_2=null", "msfragger.search_enzyme_name_2=" + protease2) \
                                         .replace("msfragger.search_enzyme_cut_2=", "msfragger.search_enzyme_cut_2=" + cuts) \
                                         .replace("msfragger.search_enzyme_nocut_2=", "msfragger.search_enzyme_nocut_2=" + nocuts) \
                                         .replace("msfragger.search_enzyme_sense_2=C", "msfragger.search_enzyme_sense_2=" + sense)
        # Specify whether to perform peptide quantification
        workflow_data = workflow_data.replace("diann.run-dia-nn=", "diann.run-dia-nn=" + str(is_perform_quantification).lower())
        workflow_filename = f"{protease1}.workflow" if not protease2 else f"{protease1}_{protease2}.workflow"
        with open(workflow_filename, 'w') as f:
            _ = f.write(workflow_data)
        logging.info(f"Workflow file generated at '{workflow_filename}'.")
        workflow_path = workflow_filename
        return True
    except Exception as error:
        logging.error("Error: An exception occurred when preparing the workflow file.")
        logging.error(error)
        return False

# Run FragPipe to perform the peptide search
def run_fragpipe(fragpipe_bin_path, num_threads, memory_limit, manifest_path, workflow_path, tools_folder_path, diann_path, python_path):
    try:
        logging.info(f"Manifest path: {manifest_path}")
        logging.info(f"Workflow path: {workflow_path}")
        logging.info("Running FragPipe...")
        fragpipe_args = [fragpipe_bin_path, "--headless", "--threads", str(num_threads), "--ram", str(memory_limit), "--manifest", manifest_path, "--workflow", workflow_path, "--workdir", os.getcwd(),
                         "--config-tools-folder", tools_folder_path, "--config-python", python_path] + ["--config-diann", diann_path] * bool(diann_path)
        run_command(fragpipe_args)
        logging.info("FragPipe finished running.")
        return True
    except Exception as error:
        logging.error("Error: An exception occurred when running FragPipe.")
        logging.error(error)
        return False

def main():
    parser = argparse.ArgumentParser()

    # required arguments
    parser.add_argument("--db_path", help = "path to the proteome database for the peptide search", type = str, required = True)
    parser.add_argument("--mass_spec_info_path", help = f"path to a mass spectrometry file list; each line is a unique mass spectrometry data file and its data type ({' / '.join(accepted_mass_spec_data_types)}, case-insensitive) separated by a tab character", type = str, required = True)
    parser.add_argument("--protease1", help = f"the first protease used", type = str, required = True, choices = protease_rules.keys())
    parser.add_argument("--output_dir", help = "output directory", type = str, required = True)
    parser.add_argument("--fragpipe_path", help = "path to the FragPipe installation directory", type = str, required = True)

    # optional arguments
    parser.add_argument("--protease2", help = f"the second protease used (optional; must not be identical to the first protease)", type = str, default = '', choices = protease_rules.keys())
    parser.add_argument("--add_contaminants", help = "add contaminants into the proteome database?", default = False, action = argparse.BooleanOptionalAction)
    parser.add_argument("--perform_quantification", help = "perform peptide quantification after the peptide search?", default = False, action = argparse.BooleanOptionalAction)

    parser.add_argument("--num_threads", help = "number of CPU threads to use for FragPipe (default: 1; specify '0' for FragPipe to use <number of cores in system - 1> threads)", default = 1, type = int)
    parser.add_argument("--memory_limit", help = "memory limit in GB for FragPipe (default: 15 GB; specify '0' to let FragPipe decide; FragPipe will exit and fail to complete the peptide search if the memory limit is below the minimum needed)", default = 15, type = int)

    parser.add_argument("--diann_path", help = "path to the DIA-NN binary for peptide quantification (default: '<fragpipe_path>/tools/diann/*/linux/diann-*')", type = str)
    parser.add_argument("--python_path", help = "path to the Python binary with easypqp installed (default: the Python binary used to run this script)", type = str)
    parser.add_argument("--philosopher_path", help = "path to the Philosopher binary for preparing the custom proteome database (default: '<fragpipe_path>/tools/Philosopher/philosopher-*')", type = str)
    parser.add_argument("--tools_folder_path", help = "path to the FragPipe tools folder containing MSFragger, IonQuant and dirTracer (default: '<fragpipe_path>/tools/')", type = str)

    args = parser.parse_args()

    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_filename = f"{timestamp}_proteomics_module.log"
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(message)s", handlers = [logging.FileHandler(log_filename), logging.StreamHandler()])
    log_file_path = get_real_path(log_filename)

    # Store the input arguments
    db_path, is_add_contaminants = args.db_path, args.add_contaminants
    mass_spec_info_path, protease1, protease2 = args.mass_spec_info_path, args.protease1.lower().strip(), args.protease2.lower().strip()
    is_perform_quantification, output_dir = args.perform_quantification, args.output_dir
    num_threads, memory_limit = args.num_threads, args.memory_limit
    fragpipe_path, diann_path, python_path, philosopher_path, tools_folder_path = args.fragpipe_path, args.diann_path, args.python_path, args.philosopher_path, args.tools_folder_path

    if not is_perform_quantification:
        diann_path = ''

    fragpipe_bin_path = os.path.join(fragpipe_path, "bin", "fragpipe")  # path to the FragPipe binary

    if not python_path:
        python_path = sys.executable

    diann_path_glob = fragpipe_path + "/tools/diann/*/linux/diann-*"
    philosopher_path_glob = fragpipe_path + "/tools/Philosopher/philosopher-*"
    diann_paths = glob.glob(diann_path_glob)
    philosopher_paths = glob.glob(philosopher_path_glob)

    if is_perform_quantification and not diann_path and len(diann_paths) != 0:
        diann_path = diann_paths[0]

    if not philosopher_path and len(philosopher_paths) != 0:
        philosopher_path = philosopher_paths[0]

    if not tools_folder_path:
        tools_folder_path = os.path.join(fragpipe_path, "tools")

    # Check the input arguments
    error_msg = ''
    if not os.path.exists(db_path):
        error_msg = f"The proteome database '{db_path}' does not exist."
    elif not os.path.isfile(db_path):
        error_msg = f"The proteome database '{db_path}' exists but it is not a file."
    elif os.path.getsize(db_path) == 0:
        error_msg = f"The proteome database '{db_path}' must not be empty."
    elif not os.path.exists(mass_spec_info_path):
        error_msg = f"The mass spectrometry file list '{mass_spec_info_path}' does not exist."
    elif not os.path.isfile(mass_spec_info_path):
        error_msg = f"The mass spectrometry file list '{mass_spec_info_path}' exists but it is not a file."
    elif os.path.getsize(mass_spec_info_path) == 0:
        error_msg = f"The mass spectrometry file list '{mass_spec_info_path}' must not be empty."
    elif protease1 not in protease_rules:
        error_msg = f"The first protease {protease1} is not supported. Please specify one from the following: {', '.join(protease_rules)}"
    elif protease2 and protease2 not in protease_rules:
        error_msg = f"The second protease {protease2} is not supported. Please specify one from the following: {', '.join(protease_rules)}"
    elif protease2 and protease1 == protease2:
        error_msg = f"Both proteases specified are identical. Please ensure both proteases are unique and come from the following: {', '.join(protease_rules)}"
    elif os.path.exists(output_dir) and not os.path.isdir(output_dir):
        error_msg = f"The output directory '{output_dir}' exists but it is not a directory."
    elif num_threads < 0:
        error_msg = f"Please specify either a positive number of CPU threads or '0' for FragPipe to use <core_number> - 1 threads."
    elif memory_limit < 0:
        error_msg = f"Please specify either a positive amount of memory to use in GB or '0' to let FragPipe decide."
    elif not os.path.exists(fragpipe_path):
        error_msg = f"The FragPipe directory '{fragpipe_path}' does not exist."
    elif not os.path.isdir(fragpipe_path):
        error_msg = f"The FragPipe directory '{fragpipe_path}' exists but it is not a directory."
    elif not os.path.exists(fragpipe_bin_path):
        error_msg = f"The FragPipe binary '{fragpipe_bin_path}' does not exist. Please ensure that the FragPipe directory specified is correct."
    elif not os.path.isfile(fragpipe_bin_path):
        error_msg = f"The FragPipe binary '{fragpipe_bin_path}' exists but it is not a file. Please ensure that the FragPipe directory specified is correct."
    elif not os.access(fragpipe_bin_path, os.X_OK):
        error_msg = f"The FragPipe binary '{fragpipe_bin_path}' is not executable. Please ensure that the FragPipe directory specified is correct and the FragPipe binary has executable permissions."
    elif os.path.getsize(fragpipe_bin_path) == 0:
        error_msg = f"The FragPipe binary '{fragpipe_bin_path}' must not be empty. Please ensure that the FragPipe directory specified is correct."
    elif is_perform_quantification and not diann_path:
        error_msg = f"There is no DIA-NN binary fitting the glob pattern '{diann_path_glob}'. Please ensure that the FragPipe directory specified is correct. Alternatively, specify the DIA-NN binary directly with --diann_path."
    elif is_perform_quantification and not os.path.exists(diann_path):
        error_msg = f"The DIA-NN binary '{diann_path}' does not exist."
    elif is_perform_quantification and not os.path.isfile(diann_path):
        error_msg = f"The DIA-NN binary '{diann_path}' exists but it is not a file."
    elif is_perform_quantification and not os.access(diann_path, os.X_OK):
        error_msg = f"The DIA-NN binary '{diann_path}' is not executable."
    elif is_perform_quantification and os.path.getsize(diann_path) == 0:
        error_msg = f"The DIA-NN binary '{diann_path}' must not be empty."
    elif not python_path:
        error_msg = f"Unable to automatically determine the location of the Python binary. Please specify the path to the Python binary directly with --python_path."
    elif not os.path.exists(python_path):
        error_msg = f"The Python binary '{python_path}' does not exist."
    elif not os.path.isfile(python_path):
        error_msg = f"The Python binary '{python_path}' exists but it is not a file."
    elif not os.access(python_path, os.X_OK):
        error_msg = f"The Python binary '{python_path}' is not executable."
    elif os.path.getsize(python_path) == 0:
        error_msg = f"The Python binary '{python_path}' must not be empty."
    elif not philosopher_path:
        error_msg = f"There is no Philosopher binary fitting the glob pattern '{philosopher_path_glob}'. Please ensure that the FragPipe directory specified is correct. Alternatively, specify the path to the Philosopher binary directly with --philosopher_path."
    elif not os.path.exists(philosopher_path):
        error_msg = f"The Philosopher binary '{philosopher_path}' does not exist."
    elif not os.path.isfile(philosopher_path):
        error_msg = f"The Philosopher binary '{philosopher_path}' exists but it is not a file."
    elif not os.access(philosopher_path, os.X_OK):
        error_msg = f"The Philosopher binary '{philosopher_path}' is not executable."
    elif os.path.getsize(philosopher_path) == 0:
        error_msg = f"The Philosopher binary '{philosopher_path}' must not be empty."
    elif not os.path.exists(tools_folder_path):
        error_msg = f"The FragPipe tools folder '{tools_folder_path}' does not exist."
        if not args.tools_folder_path:
            error_msg += " Please ensure that the FragPipe directory specified is correct. Alternatively, specify the path to the FragPipe tools folder directly with --tools_folder_path."
    elif not os.path.isdir(tools_folder_path):
        error_msg = f"The FragPipe tools folder '{tools_folder_path}' exists but it is not a directory."
        if not args.tools_folder_path:
            error_msg += " Please ensure that the FragPipe directory specified is correct. Alternatively, specify the path to the FragPipe tools folder directly with --tools_folder_path."

    if error_msg:
        logging.error(f"Error: {error_msg}")
        sys.exit(1)

    # Turn all paths provided into real paths
    db_path = get_real_path(db_path)
    output_dir = get_real_path(output_dir)
    fragpipe_bin_path = get_real_path(fragpipe_bin_path)
    if diann_path:
        diann_path = get_real_path(diann_path)
    python_path = get_real_path(python_path)
    philosopher_path = get_real_path(philosopher_path)
    tools_folder_path = get_real_path(tools_folder_path)

    # Parse the mass spectrometry file list
    mass_spec_filenames = []
    mass_spec_data_types = []
    try:
        error_msg = ''
        with open(mass_spec_info_path, 'r') as f:
            for line in f:
                cols = [col.strip() for col in line.strip().split('\t')]
                if len(cols) != 2:
                    error_msg = "Each line of the mass spectrometry file list must consist of two columns (the mass spectrometry data file and its data type) separated by a tab character."
                    break
                [mass_spec_filename, mass_spec_data_type] = cols
                if not os.path.exists(mass_spec_filename):
                    error_msg = f"The mass spectrometry data file '{mass_spec_filename}' does not exist. Please ensure all data files specified in the mass spectrometry file list are correct."
                    break
                if not os.path.isfile(mass_spec_filename):
                    error_msg = f"The mass spectrometry data file '{mass_spec_filename}' exists but it is not a file. Please ensure all data files specified in the mass spectrometry file list are correct."
                    break
                if os.path.getsize(mass_spec_filename) == 0:
                    error_msg = f"The mass spectrometry data file '{mass_spec_filename}' is empty. Please ensure all data files specified in the mass spectrometry file list are not empty."
                if mass_spec_data_type.lower() not in accepted_mass_spec_data_types_lower:
                    error_msg = f"The mass spectrometry data type '{mass_spec_data_type}' for the file '{mass_spec_filename}' is not supported. Please ensure all data types specified in the mass spectrometry file list are within the following group (case-insensitive): {', '.join(accepted_mass_spec_data_types)}"
                    break
                # Enforce the casing FragPipe uses for the mass spectrometry data types
                for accepted_mass_spec_data_type in accepted_mass_spec_data_types:
                    if mass_spec_data_type.lower() == accepted_mass_spec_data_type.lower():
                        mass_spec_data_type = accepted_mass_spec_data_type
                        break
                mass_spec_filenames.append(mass_spec_filename)
                mass_spec_data_types.append(mass_spec_data_type)
        if len(mass_spec_filenames) != len(set(mass_spec_filenames)):
            error_msg = "Some data files specified in the mass spectrometry file list appeared more than once. Please ensure all data files specified are unique."
        if error_msg:
            logging.error(f"Error: {error_msg}")
            sys.exit(1)
    except Exception as error:
        logging.error("Error: An exception occurred when parsing the mass spectrometry file list.")
        logging.error(error)
        sys.exit(1)

    # Print the user configuration options
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Proteome database path: {db_path}")
    if protease2:
        logging.info(f"Proteases: {protease1}, {protease2}")
    else:
        logging.info(f"Protease: {protease1}")
    logging.info(f"Add contaminants: {'Yes' if is_add_contaminants else 'No'}")
    logging.info(f"Perform peptide quantification: {'Yes' if is_perform_quantification else 'No'}")
    logging.info(f"CPU threads: {'<number of cores in system - 1>' if num_threads == 0 else num_threads}")
    logging.info(f"Memory limit: {'Let FragPipe decide' if memory_limit == 0 else str(memory_limit) + ' GB'}")
    logging.info(f"Mass spectrometry data files: {', '.join(mass_spec_filenames)}")
    logging.info(f"Mass spectrometry data types: {', '.join(mass_spec_data_types)}")
    logging.info(f"FragPipe binary path: {fragpipe_bin_path}")
    logging.info(f"FragPipe tools folder path: {tools_folder_path}")
    logging.info(f"Philosopher binary path: {philosopher_path}")
    logging.info(f"Python binary path: {python_path}")
    if is_perform_quantification:
        logging.info(f"DIA-NN binary path: {diann_path}")

    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created the output directory '{output_dir}'.")

    # Go into the output directory
    os.chdir(output_dir)

    # Prepare the input files and run FragPipe
    if prepare_proteome_database(db_path, philosopher_path, is_add_contaminants) and prepare_manifest(mass_spec_filenames, mass_spec_data_types) and prepare_workflow(protease1, protease2, is_perform_quantification) and run_fragpipe(fragpipe_bin_path, num_threads, memory_limit, manifest_path, workflow_path, tools_folder_path, diann_path, python_path):
        # Find the FragPipe log file
        fragpipe_log_files = glob.glob(os.path.join(output_dir, "log_*.txt"))
        if len(fragpipe_log_files) > 0:
            fragpipe_log_file = fragpipe_log_files[0]
            logging.info(f"The complete FragPipe run log is located at {fragpipe_log_file}.")

        # Find the peptide search results
        is_error = False
        fragpipe_peptide_tsv = os.path.join(output_dir, "peptide.tsv")
        if os.path.exists(fragpipe_peptide_tsv):
            logging.info(f"The peptide search results are located at {fragpipe_peptide_tsv}.")
        else:
            logging.error(f"Error: Peptide search results are not located at {fragpipe_peptide_tsv}.")
            is_error = True

        # Find the peptide quantification results
        if is_perform_quantification:
            fragpipe_report_pr_matrix_tsv = os.path.join(output_dir, "dia-quant-output", "report.pr_matrix.tsv")
            if os.path.exists(fragpipe_report_pr_matrix_tsv):
                logging.info(f"The peptide quantification results are located at {fragpipe_report_pr_matrix_tsv}.")
            else:
                logging.error(f"Error: Peptide quantification results are not located at {fragpipe_report_pr_matrix_tsv}.")
                is_error = True

        if not is_error:
            logging.info("Proteomics search complete.")
        else:
            if len(fragpipe_log_files) > 0:
                logging.error(f"Proteomics search failed. Please read the log file at {log_file_path} and the complete FragPipe run log for details on the errors encountered.")
            else:
                logging.error(f"Proteomics search failed. Please read the log file at {log_file_path} for details on the errors encountered.")
            sys.exit(1)
    else:
        logging.error(f"Proteomics search failed. Please read the log file at {log_file_path} for details on the errors encountered.")
        sys.exit(1)

if __name__ == "__main__":
    main()
