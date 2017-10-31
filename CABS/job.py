"""
Module for running cabsDock jobs.
"""
import operator
import os

from abc import ABCMeta, abstractmethod
from CABS import logger, PDBlib, cabs
from CABS.align import save_csv
from CABS.cluster import Clustering
from CABS.cmap import ContactMapFactory
from CABS.filter import Filter
from CABS.plots import graph_RMSF, plot_E_RMSD, plot_RMSD_N
from CABS.protein import ProteinComplex
from CABS.restraints import Restraints
from CABS.trajectory import Trajectory
from CABS.utils import SCModeler, CONFIG_HEADER


_name = 'JOB'


class CABSTask(object):
    """Abstract CABS job instance."""

    __metaclass__ = ABCMeta

    def __init__(self, **kwargs):

        # self.__dict__.update(kwargs)
        self.aa_rebuild = kwargs.get('aa_rebuild')
        self.add_peptide = kwargs.get('add_peptide')
        self.align = kwargs.get('align')
        self.align_options = dict(kwargs.get('align_options', []))
        self.align_peptide_options = dict(kwargs.get('align_peptide_options', []))
        self.ca_rest_add = kwargs.get('ca_rest_add')
        self.ca_rest_file = kwargs.get('ca_rest_file')
        self.ca_rest_weight = kwargs.get('ca_rest_weight')
        self.clustering_iterations = kwargs.get('clustering_iterations')
        self.clustering_medoids = kwargs.get('clustering_medoids')
        self.contact_maps = kwargs.get('contact_maps')
        self.contact_threshold = kwargs.get('contact_threshold')
        self.dssp_command = kwargs.get('dssp_command')
        self.exclude = kwargs.get('exclude')
        self.excluding_distance = kwargs.get('excluding_distance')
        self.filtering_count = kwargs.get('filtering_count')
        self.filtering_mode = kwargs.get('filtering_mode')
        self.fortran_command = kwargs.get('fortran_command')
        self.image_file_format = kwargs.get('image_file_format')
        self.input_protein = kwargs.get('input_protein')
        self.insertion_attempts = kwargs.get('insertion_attempts')
        self.insertion_clash = kwargs.get('insertion_clash')
        self.load_cabs_files = kwargs.get('load_cabs_files')
        self.mc_annealing = kwargs.get('mc_annealing')
        self.mc_cycles = kwargs.get('mc_cycles')
        self.mc_steps = kwargs.get('mc_steps')
        self.modeller_iterations = kwargs.get('modeller_iterations')
        self.pdb_output = kwargs.get('pdb_output')
        self.peptide = kwargs.get('peptide')
        self.protein_flexibility = kwargs.get('protein_flexibility')
        self.protein_restraints = kwargs.get('protein_restraints')
        self.protein_restraints_reduce = kwargs.get('protein_restraints_reduce')
        self.random_seed = kwargs.get('random_seed')
        self.reference_pdb = kwargs.get('reference_pdb')
        self.replicas = kwargs.get('replicas')
        self.replicas_dtemp = kwargs.get('replicas_dtemp')
        self.save_cabs_files = kwargs.get('save_cabs_files')
        self.save_config = kwargs.get('save_config')
        self.sc_rest_add = kwargs.get('sc_rest_add')
        self.sc_rest_file = kwargs.get('sc_rest_file')
        self.sc_rest_weight = kwargs.get('sc_rest_weight')
        self.separation = kwargs.get('separation')
        self.temperature = kwargs.get('temperature')
        self.verbose = kwargs.get('verbose')
        self.work_dir = kwargs.get('work_dir')
        self.weighted_fit = kwargs.get('weighted_fit')
        self.weighted_fit_file = kwargs.get('weighted_fit_file')
        self.colors = kwargs.get('contact_map_colors')

        # Job attributes collected.
        self.config = kwargs
        self.initial_complex = None
        self.restraints = None
        self.cabsrun = None
        self.trajectory = None
        self.filtered_trajectory = None
        self.filtered_ndx = None
        self.medoids = None
        self.clusters_dict = None
        self.clusters = None
        self.rmslst = {}
        self.results = None
        self.reference = None

        # Workdir processing: making sure work_dir is abspath
        self.work_dir = os.path.abspath(self.work_dir)
        logger.setup_log_level(self.verbose)
        try:
            os.makedirs(self.work_dir)
        except OSError:
            if os.path.isdir(self.work_dir):
                logger.warning(_name, '{} already exists. Output data will be overwritten.'.format(self.work_dir))
            else:
                logger.exit_program(
                    _name, '{} already exists and is not a directory. Choose different name.'.format(self.work_dir)
                )

        try:
            self.file_TRAF, self.file_SEQ = self.load_cabs_files
        except (ValueError, TypeError):
            self.file_TRAF = self.file_SEQ = None

        # self.peptide + self.add_peptide -> self.ligand
        self.peptides = []
        if self.peptide:
            self.peptides.extend([[p, 'random', 'random'] for p in self.peptide])
        if self.add_peptide:
            self.peptides.extend([p for p in self.add_peptide if p])

        # Pdb output processing
        if 'A' in self.pdb_output:
            self.pdb_output = 'RFCM'
        elif 'N' in self.pdb_output:
            self.pdb_output = ''

        if self.weighted_fit_file:
            self.weighted_fit = self.weighted_fit_file

        _DSSP_COMMAND = self.dssp_command
        _FORTRAN_COMMAND = self.fortran_command

    def run(self):
        ftraf = self.file_TRAF
        fseq = self.file_SEQ
        self.setup_job()
        if self.reference_pdb:
            self.parse_reference(self.reference_pdb)
        withcabs = None in (ftraf, fseq)

        if withcabs:
            self.setup_cabs_run()
            self.execute_cabs_run()
        self.load_output(ftraf, fseq)
        self.score_results(
            n_filtered=self.filtering_count,
            number_of_medoids=self.clustering_medoids,
            number_of_iterations=self.clustering_iterations
        )

        if self.reference_pdb:
            self.calculate_rmsd()
        self.save_config_file()
        self.draw_plots(colors=self.colors)
        self.save_models()
        logger.info(module_name=_name, msg='Simulation completed successfully')

    @abstractmethod
    def setup_job(self):
        pass

    @abstractmethod
    def calculate_rmsd(self):
        pass

    @abstractmethod
    def parse_reference(self, ref):
        pass

    @abstractmethod
    def draw_plots(self, plots_dir=None, colors=['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000']):
        pass

    @abstractmethod
    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, plots_dir):
        scmodeler = SCModeler(ca_traj.template)
        sc_traj_full = scmodeler.calculate_sc_traj(ca_traj.coordinates)
        sc_med = scmodeler.calculate_sc_traj(meds.coordinates)
        shp = sc_med.shape
        sc_med = sc_med.reshape((shp[1], shp[0]) + shp[2:])

        cmapdir = os.path.join(self.work_dir, 'contact_maps')
        try:
            os.mkdir(cmapdir)
        except OSError:
            pass

        return sc_traj_full, sc_med, cmapdir

    def prepare_restraints(self):

        # generate protein restraints
        protein_restraints = Restraints(
            self.initial_complex.protein.generate_restraints(*self.protein_restraints)
        )

        # reduce number of restraints
        if self.protein_restraints_reduce:
            protein_restraints.reduce_by(self.protein_restraints_reduce)

        # additional restraints
        add_restraints = Restraints('')

        if self.ca_rest_add:
            add_restraints += Restraints.from_parser(self.ca_rest_add)

        if self.sc_rest_add:
            add_restraints += Restraints.from_parser(self.sc_rest_add, sg=True)

        if self.ca_rest_file:
            for filename in self.ca_rest_file:
                add_restraints += Restraints.from_file(filename)

        if self.sc_rest_file:
            for filename in self.sc_rest_file:
                add_restraints += Restraints.from_file(filename, sg=True)

        protein_restraints += add_restraints.update_id(self.initial_complex.new_ids)
        return protein_restraints

    def save_config_file(self):
        if self.save_config:
            with open(os.path.join(self.work_dir, 'config.ini'), 'w') as configfile:
                configfile.write(CONFIG_HEADER)
                for k in sorted(self.config):
                    peptide_counter = 0
                    value = self.config[k]
                    line = str(k)+': '
                    if k == 'ligand':
                        for lgnd in value:
                            if peptide_counter == 0:
                                line = 'peptide: '
                            else:
                                line += '\nadd-peptide: '
                            for element in lgnd:
                                line += str(element)+' '
                            peptide_counter += 1
                    elif isinstance(value, tuple):
                        line = str(k) + ': '
                        for item in value:
                            line += str(item)+' '
                    else:
                        line = str(k) + ': '+str(value)
                    configfile.write('\n'+line)

    def setup_cabs_run(self):
        logger.info(module_name="CABS", msg='Setting up CABS simulation.')
        # Initializing CabsRun instance
        self.cabsrun = cabs.CabsRun(
            protein_complex=self.initial_complex,
            restraints=self.prepare_restraints(),
            work_dir=self.work_dir,
            excluding_distance=self.excluding_distance,
            replicas=self.replicas,
            replicas_dtemp=self.replicas_dtemp,
            temperature=self.temperature,
            ca_rest_weight=self.ca_rest_weight,
            sc_rest_weight=self.sc_rest_weight,
            mc_annealing=self.mc_annealing,
            mc_cycles=self.mc_cycles,
            mc_steps=self.mc_steps
        )
        return self.cabsrun

    def execute_cabs_run(self):
        self.cabsrun.run()

    @abstractmethod
    def load_output(self, ftraf=None, fseq=None):
        """
        Method for loading previously done simulation results. Stores the results to self.trajectory.
        :param ftraf: path to TRAF file
        :param fseq: path to SEQ file
        :return: returns trajectory.Trajectory instance
        """
        if ftraf is not None and fseq is not None:
            logger.debug(module_name=_name, msg="Loading trajectories from: %s, %s" % (ftraf, fseq))
            self.trajectory = Trajectory.read_trajectory(ftraf, fseq)
        else:
            logger.debug(module_name=_name, msg="Loading trajectories from the CABS run")
            self.trajectory = self.cabsrun.get_trajectory()
        self.trajectory.weights = self.initial_complex.protein.weights
        self.trajectory.template.update_ids(self.initial_complex.protein.old_ids, pedantic=False)
        chs = ''.join(self.initial_complex.protein_chains)
        tchs = ''.join(set(chs).intersection(self.trajectory.template.list_chains().keys()))
        self.trajectory.tmp_target_chs = tchs
        # TODO: target chains of trajectory template... should not be added here.
        # or even: should not be needed to add at all.
        ic_stc, tt_stc, dummy_aln = self.trajectory.align_to(self.initial_complex.protein, chs, tchs, align_mth='SW')
        self.trajectory.superimpose_to(ic_stc, tt_stc)
        logger.info(module_name=_name, msg="Trajectories loaded successfully")
        return self.trajectory

    @abstractmethod
    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        pass

    def save_models(self):
        # output folder
        output_folder = os.path.join(self.work_dir, 'output_pdbs')
        logger.log_file(module_name=_name, msg="Saving pdb files to " + str(output_folder))
        try:
            os.mkdir(output_folder)
        except OSError:
            logger.warning(_name, "Possibly overwriting previous pdb files. Use --work-dir <DIR> to avoid that.")
        # Saving the trajectory to PDBs:
        if 'R' in self.pdb_output:
            logger.log_file(module_name=_name, msg='Saving replicas...')
            self.trajectory.to_pdb(mode='replicas', to_dir=output_folder)
        # Saving top1000 models to PDB:
        if 'F' in self.pdb_output:
            logger.log_file(module_name=_name, msg='Saving filtered models...')
            self.filtered_trajectory.to_pdb(mode='replicas', to_dir=output_folder, name='top1000')
        # Saving clusters in CA representation
        if 'C' in self.pdb_output:
            logger.log_file(module_name=_name, msg='Saving clusters...')
            for i, cluster in enumerate(self.clusters):
                cluster.to_pdb(mode='replicas', to_dir=output_folder, name='cluster_{0}'.format(i))
        # Saving final models:
        if 'M' in self.pdb_output:
            if self.aa_rebuild:
                logger.log_file(module_name=_name, msg='Saving final models (in AA representation)')
                pdb_medoids = self.medoids.to_pdb()
                from CABS.ca2all import ca2all
                for i, fname in enumerate(pdb_medoids):
                    ca2all(
                        fname,
                        output=os.path.join(output_folder, 'model_{0}.pdb'.format(i)),
                        iterations=self.modeller_iterations,
                        out_mdl=os.path.join(self.work_dir, 'output_data','modeller_output_{0}.txt'.format(i)),
                        work_dir=self.work_dir
                    )
            else:
                logger.log_file(module_name=_name, msg='Saving final models (in CA representation)')
                self.medoids.to_pdb(mode='models', to_dir=output_folder, name='model')


class DockTask(CABSTask):
    """Class representing single CABS job."""

    def setup_job(self):
        if not self.peptides:
            raise ValueError('No peptide given')
        self.initial_complex = ProteinComplex(
            protein=self.input_protein,
            flexibility=self.protein_flexibility,
            exclude=self.exclude,
            weights=self.weighted_fit,
            peptides=self.peptides,
            replicas=self.replicas,
            separation=self.separation,
            insertion_attempts=self.insertion_attempts,
            insertion_clash=self.insertion_clash,
            work_dir=self.work_dir,
            dssp_command=self.dssp_command
        )

    def load_output(self, ftraf=None, fseq=None):
        """
        Method for loading previously done simulation results. Stores the results to self.trajectory.
        :param ftraf: path to TRAF file
        :param fseq: path to SEQ file
        :return: returns trajectory.Trajectory instance
        """
        ret = super(DockTask, self).load_output(ftraf, fseq)
        ret.number_of_peptides = len(self.peptides)
        return ret

    def calculate_rmsd(self, save=True):
        logger.debug(module_name=_name, msg="RMSD calculations starting...")
        if save:
            odir = os.path.join(self.work_dir, 'output_data')
            try:
                os.mkdir(odir)
            except OSError:
                pass
        all_results = {}
        ref_trg_stc, self_trg_stc, trg_aln = self.trajectory.align_to(
            self.reference[0], self.reference[1], self.trajectory.tmp_target_chs,
            align_mth=self.align, kwargs=self.align_options
        )
        self.trajectory.superimpose_to(ref_trg_stc, self_trg_stc)
        if save:
            sfname = os.path.join(self.work_dir, 'output_data', 'reference_alignment')
            paln_trg = sfname + '_target.csv'
            save_csv(paln_trg, ('reference', 'template'), trg_aln)
        for pept_chain, ref_pept_chain in zip(self.initial_complex.peptide_chains, self.reference[2]):
            ref_pep_stc, self_pep_stc, pep_aln = self.trajectory.align_to(
                self.reference[0], ref_pept_chain, pept_chain, align_mth=self.align,
                kwargs=self.align_peptide_options
            )
            if save:
                paln_pep = sfname + '_%s.csv' % pept_chain
                save_csv(paln_pep, ('reference', 'template'), pep_aln)
            self.rmslst[pept_chain] = self.trajectory.rmsd_to_reference(ref_pep_stc, self_pep_stc)
            rmsds = [header.rmsd for header in self.medoids.headers]
            results = {}
            results['rmsds_all'] = [header.rmsd for header in self.trajectory.headers]
            results['rmsds_filtered'] = [header.rmsd for header in self.filtered_trajectory.headers]
            results['rmsds_medoids'] = rmsds
            results['lowest_all'] = sorted(results['rmsds_all'])[0]
            results['lowest_filtered'] = sorted(results['rmsds_filtered'])[0]
            results['lowest_medoids'] = sorted(results['rmsds_medoids'])[0]
            # Saving rmsd results
            if save:
                with open(os.path.join(odir, 'lowest_rmsds_%s.txt' % pept_chain), 'w') as outfile:
                    outfile.write(
                        'lowest_all; lowest_filtered; lowest_medoids\n {0};{1};{2}'.format(
                            results['lowest_all'], results['lowest_filtered'], results['lowest_medoids']
                        )
                    )
                for _type in ['all', 'filtered', 'medoids']:
                    with open(os.path.join(odir, '{0}_rmsds_{1}.txt'.format(_type, pept_chain)), 'w') as outfile:
                        for rmsd in results['rmsds_' + _type]:
                            outfile.write(str(rmsd) + ';\n')
            all_results[pept_chain] = results
        logger.info(module_name=_name, msg="RMSD successfully saved")
        return all_results

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        logger.debug(module_name=_name, msg="Scoring results")
        # Filtering the trajectory
        self.filtered_trajectory, self.filtered_ndx = Filter(self.trajectory, n_filtered).cabs_filter()
        # Clustering the trajectory
        self.medoids, self.clusters_dict, self.clusters = Clustering(
            self.filtered_trajectory,
            'chain ' + ','.join(
                self.initial_complex.peptide_chains,
            )
        ).cabs_clustering(number_of_medoids=number_of_medoids, number_of_iterations=number_of_iterations)
        logger.info(module_name=_name, msg="Scoring results successful")

    def draw_plots(self, plots_dir=None, colors=['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000']):
        logger.debug(module_name=_name, msg="Drawing plots")
        # set the plots dir
        if plots_dir is None:
            pltdir = os.path.join(self.work_dir, 'plots')
            try:
                os.mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir
        logger.log_file(module_name=_name, msg="Saving plots to %s" % pltdir)

        graph_RMSF(self.trajectory, self.initial_complex.protein_chains, os.path.join(pltdir, 'RMSF'))

        # RMSD-based graphs
        if self.reference_pdb:
            logger.log_file(module_name=_name, msg="Saving RMSD plots")
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD(
                    [self.trajectory, self.filtered_trajectory],
                    [rmslst, rmslst[self.filtered_ndx, ]],
                    ['all models', 'top 1000 models'],
                    os.path.join(pltdir, 'E_RMSD_%s' % k)
                )
                plot_RMSD_N(
                    rmslst.reshape(self.replicas, -1),
                    os.path.join(pltdir, 'RMSD_frame_%s' % k)
                )

        # Contact maps
        if self.contact_maps:
            logger.log_file(module_name=_name, msg="Saving contact maps")
            self.mk_cmaps(
                self.trajectory, self.medoids, self.clusters_dict, self.filtered_ndx, 4.5, pltdir, colors=colors
            )
        logger.info(module_name=_name, msg="Plots successfully saved")

    def mk_cmaps(
            self, ca_traj, meds, clusts, top1k_inds, thr, plots_dir,
            colors=['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000']
    ):
        sc_traj_full, sc_med, cmapdir = super(DockTask, self).mk_cmaps(
            ca_traj, meds, clusts, top1k_inds, thr, plots_dir
        )

        sc_traj_1k = sc_traj_full.reshape(1, -1, len(ca_traj.template), 3)[:, top1k_inds, :, :]

        rchs = self.initial_complex.protein_chains
        lchs = self.initial_complex.peptide_chains

        targ_cmf = ContactMapFactory(rchs, rchs, ca_traj.template)
        cmfs = {lig: ContactMapFactory(rchs, lig, ca_traj.template) for lig in lchs}
        cmap10ktarg = reduce(operator.add, targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg.zero_diagonal()
        cmap10ktarg.save_all(cmapdir + '/target_all', break_long_x=0, norm_n=True, colors=colors)

        for lig, cmf in cmfs.items():
            cmaps = cmf.mk_cmap(sc_traj_full, thr)
            for n, cmap in enumerate(cmaps):
                cmap.save_all(cmapdir + '/replica_%i_ch_%s' % (n + 1, lig), norm_n=True, colors=colors)
            cmap10k = reduce(operator.add, cmaps)
            cmap10k.save_all(cmapdir + '/all_ch_%s' % lig, norm_n=True, colors=colors)
            cmap10k.save_histo(plots_dir + '/all_contacts_histo_%s' % lig)
            cmap1k = cmf.mk_cmap(sc_traj_1k, thr)[0]
            cmap1k.save_all(cmapdir + '/top1000_ch_%s' % lig, norm_n=True, colors=colors)
            cmaps_top = cmf.mk_cmap(sc_med, thr)
            for n, cmap in enumerate(cmaps_top):
                cmap.save_all(cmapdir + '/top_%i_ch_%s' % (n + 1, lig), norm_n=True, colors=colors)
            for cn, clust in clusts.items():
                ccmap = cmf.mk_cmap(sc_traj_1k, thr, frames=clust)[0]
                ccmap.save_all(cmapdir + '/cluster_%i_ch_%s' % (cn, lig), norm_n=True, colors=colors)

    def parse_reference(self, ref):
        try:
            source, rec, pep = ref.split(':')
            self.reference = (PDBlib.Pdb(ref, selection='name CA', no_exit=True, verify=True).atoms, rec, pep)
            if len(self.initial_complex.peptide_chains) != len(self.reference[2]):
                raise ValueError
        except (ValueError, PDBlib.Pdb.InvalidPdbInput):
            logger.warning(_name, 'Invalid reference {}'.format(ref))
            self.reference = None


class FlexTask(CABSTask):
    """Class of CABSFlex jobs."""

    def setup_job(self):
        self.initial_complex = ProteinComplex(
            protein=self.input_protein,
            flexibility=self.protein_flexibility,
            exclude=self.exclude,
            weights=self.weighted_fit,
            peptides=self.peptides,
            replicas=self.replicas,
            separation=self.separation,
            insertion_attempts=self.insertion_attempts,
            insertion_clash=self.insertion_clash,
            work_dir=self.work_dir,
            dssp_command=self.dssp_command
        )

        if self.reference_pdb is None:
            self.reference_pdb = True

        # remove filtered trajectory from pdb saving
        self.pdb_output = self.pdb_output.replace('F', '')

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        # Clustering the trajectory
        clst = Clustering(self.trajectory, 'chain ' + ','.join(self.initial_complex.protein_chains))
        self.medoids, self.clusters_dict, self.clusters = clst.cabs_clustering(number_of_medoids=number_of_medoids, number_of_iterations=number_of_iterations)
        self.rmslst = {self.initial_complex.protein_chains: clst.distance_matrix[0]}

    def load_output(self, *args, **kwargs):
        ret = super(FlexTask, self).load_output(*args, **kwargs)
        ret.number_of_peptides = 0
        return ret

    def calculate_rmsd(self, reference_pdb=None, save=True):
        
        logger.debug(module_name=_name, msg="RMSD calculations starting...")
        if save:
            odir = os.path.join(self.work_dir, 'output_data')
            try:
                os.mkdir(odir)
            except OSError:
                pass

        chs_ids = self.trajectory.tmp_target_chs

        ref_trg_stc, self_trg_stc, trg_aln = self.trajectory.align_to(
            self.reference[0], self.reference[1], chs_ids,
            align_mth=self.align, kwargs=self.align_options
        )
        self.trajectory.superimpose_to(ref_trg_stc, self_trg_stc)
        if save:
            sfname = os.path.join(self.work_dir, 'output_data', 'reference_alignment')
            paln_trg = sfname + '_target.csv'
            save_csv(paln_trg, ('reference', 'template'), trg_aln)
        self.rmslst[chs_ids] = self.trajectory.rmsd_to_reference(ref_trg_stc, self_trg_stc)
        rmsds = [header.rmsd for header in self.medoids.headers]
        results = {}
        results['rmsds_all'] = [header.rmsd for header in self.trajectory.headers]
        results['rmsds_medoids'] = rmsds
        results['lowest_all'] = sorted(results['rmsds_all'])[0]
        results['lowest_medoids'] = sorted(results['rmsds_medoids'])[0]
        # Saving rmsd results
        if save:
            with open(os.path.join(odir, 'lowest_rmsds_%s.txt' % chs_ids), 'w') as outfile:
                outfile.write(
                    'lowest_all; lowest_medoids\n {0};{1}'.format(
                        results['lowest_all'], results['lowest_medoids']
                    )
                )
            for _type in ['all', 'medoids']:
                with open(os.path.join(odir, '{0}_rmsds_{1}.txt'.format(_type, chs_ids)), 'w') as outfile:
                    for rmsd in results['rmsds_' + _type]:
                        outfile.write(str(rmsd) + ';\n')
        logger.info(module_name=_name, msg="RMSD successfully saved")
        return {chs_ids: results}

    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, plots_dir, colors=['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000']):
        sc_traj_full, sc_med, cmapdir = super(FlexTask, self).mk_cmaps(
            ca_traj, meds, clusts, top1k_inds, thr, plots_dir
        )

        rchs = self.initial_complex.protein_chains

        cmf = ContactMapFactory(rchs, rchs, ca_traj.template)
        cmap_all = reduce(operator.add, cmf.mk_cmap(sc_traj_full, thr))

        cmaptop = reduce(operator.add, cmf.mk_cmap(sc_med, thr))

        for cmap, fname in zip((cmap_all, cmaptop), ('all', 'top10')):
            cmap.zero_diagonal()
            cmap.save_all(cmapdir + '/' + fname, break_long_x=0, norm_n=True, colors=colors)

    def parse_reference(self, ref):
        try:
            try:
                dummy, trg_chids = ref.split(":")
                self.reference = (
                    PDBlib.Pdb(ref, selection='name CA', no_exit=True, verify=True).atoms, trg_chids
                )
            except AttributeError: # if ref is None it has no split mth
                self.reference = (self.initial_complex, self.initial_complex.protein_chains)
        except PDBlib.Pdb.InvalidPdbInput:
            logger.warning(_name, 'Invalid reference {}'.format(ref))

    def draw_plots(self, plots_dir=None, colors=['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000']):
        # set the plots dir
        if plots_dir is None:
            pltdir = os.path.join(self.work_dir, 'plots')
            try:
                os.mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir

        graph_RMSF(self.trajectory,
                   self.initial_complex.protein_chains,
                   os.path.join(pltdir, 'RMSF'),
                   False,
                   self.image_file_format)

        # RMSD-based graphs
        if self.reference_pdb:
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD(
                    [self.trajectory],
                    [rmslst],
                    ['all models'],
                    os.path.join(pltdir, 'E_RMSD_%s' % k),
                    self.image_file_format,
                    interaction=False
                )
                plot_RMSD_N(
                    rmslst.reshape(self.replicas, -1),
                    os.path.join(pltdir, 'RMSD_frame_%s' % k),
                    self.image_file_format
                )

        # Contact maps
        if self.contact_maps:
            self.mk_cmaps(
                self.trajectory,
                self.medoids,
                self.clusters_dict,
                self.filtered_ndx,
                self.contact_threshold,
                pltdir,
                colors=colors
            )
