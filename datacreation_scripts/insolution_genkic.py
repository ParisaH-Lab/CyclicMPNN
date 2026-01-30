#!/usr/bin/env python
#
# Author: Andrew Powers (apowers4@uoregon.edu | apowers@flatironinstitute.org)
# @brief: Script to generate cyclic peptide ensembles in solution using GenKIC.
#       For more information see https://github.com/PowersPope/UsefulRosettaScripts
#

# Packages
import argparse
from typing import List, Tuple
import datetime

# decorators
from decorators import timeit

# Rosetta imports
from pyrosetta import init
import pyrosetta.io as io
import pyrosetta.rosetta.core as core
import pyrosetta.rosetta.protocols as protocols
import pyrosetta.rosetta.protocols.generalized_kinematic_closure as genkic
import pyrosetta.rosetta.protocols.simple_moves as sm
import pyrosetta.rosetta.core.select.residue_selector as select
from pyrosetta.rosetta.protocols.cyclic_peptide import PeptideInternalHbondsFilter
import pyrosetta.rosetta.core.io.silent as silent

# Functions
def foldtree_define(N: int) -> int:
    """
    Define what the root of our fold tree is

    PARAMS
    ------
    :N (int): Desired size of input structure

    RETURNS
    -------
    :root (int): The root residue index (zero indexed)
    """
    if N % 2 == 0:
            return int(N/2)
    else:
        return int((N-1)/2)

def peptide_stub_mover(
        pose: core.pose.Pose, 
        residue_type: str,
        sequence_len: int,
        ) -> core.pose.Pose:
    """Add residues to a pose, if you don't use pose_from_sequence (You should use pose_from_sequence)

    PARAMS
    ------
    :pose: Our empty pose
    :residue_type: The residue type yo would like to add
    :sequence_len: The number of repeats that residue_type should be added

    RETURNS
    -------
    An inplace pose with residue_type added sequence_len times
    """
    # init our Mover using XML, since there are extra options that arent needed in XML
    # whereas PyRosetta wants them specified and they dont have the defaults.
    xml = f"""
    <ROSETTASCRIPTS>
        <MOVERS>
            j
            <PeptideStubMover name="extension" reset="false">
                <Append resname="{residue_type}" repeat="{sequence_len}"/>
            </PeptideStubMover>
        </MOVERS>
        <PROTOCOLS>
            <Add mover="extension"/>
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """
    pepStubMover = protocols.rosetta_scripts.XmlObjects.create_from_string(xml).get_mover("extension")
    pepStubMover.apply(pose)
    return pose.clone()

def flip_trans_omegabonds(
        pose: core.pose.Pose,
        sequence_len: int,
        ) -> int:
    """Flip the omega angles in the backbone to trans orientation
    this needs to be done if PeptideStubMover was used first

    PARAMS
    ------
    :pose: Our filled in pose
    :sequence_len: The size of our peptide

    RETURNS
    -------
    A pose with trans omega bonds, as opposed to the default cis.
    """
    # First get our root
    root = foldtree_define(sequence_len)
    # Unfortunately SetTorsion mover (which is used in RosettaScripts XML)
    # Is not fully converted well to PyRosetta, so we will specify it in XML
    xmlSetTorsion = f"""
    <ROSETTASCRIPTS>
        <MOVERS>
            <SetTorsion name="flip_trans" foldtree_root="{root}">
                <Torsion residue="{','.join([str(i) for i in range(1, sequence_len+1)])}" torsion_name="omega" angle="180"/>
            </SetTorsion>
        </MOVERS>
        <PROTOCOLS>
            <Add mover="flip_trans"/>
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """
    setTorsion = protocols.rosetta_scripts.XmlObjects.create_from_string(xmlSetTorsion).get_mover("flip_trans")
    setTorsion.apply(pose)

    return 0

# Class Generation
class BackboneGeneration:
    def __init__(
        self,
        debug: bool = False,
        randomize_root: bool = True,
        time_test: bool = False,
        empty_pose_test: bool = False,
        ):
        """
        Setup and init our Pyrosetta instance to generate macrocycle backbones + Design

        PARAMS
        ------
        :debug: Output Debug Print Messages
        :randomize_root: This should be set if your root residue (which is an anchor) needs
            to be randomized. If you had a target, you would not randomize this
        :time_test: Comparison for XML vs PyRosetta. Since XML writes to disk every struct
        :empty_pose_test: Instead of using pose_from_sequence, we will create an empty
            pose and then append to it.
        """
        # Set up our instance of PyRosetta
        if debug:
            init("-score:weights ref2015 -out:file:silent_struct_type binary -overwrite")
        else:
            init("-mute all -score:weights ref2015 -out:file:silent_struct_type binary -overwrite")

        # init class global variables
        self.scorefxn = core.scoring.get_score_function(is_fullatom = True)
        self.scorefxn.set_weight(core.scoring.hbond_lr_bb, 1.0)
        self.scorefxn.set_weight(core.scoring.hbond_sr_bb, 1.0)
        self.scorefxn.set_weight(core.scoring.coordinate_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.atom_pair_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.dihedral_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.angle_constraint, 1.0)
        self.scorefxn.set_weight(core.scoring.chainbreak, 1.0)
        emopts = core.scoring.methods.EnergyMethodOptions( self.scorefxn.energy_method_options() ) # make a copy
        emopts.hbond_options().decompose_bb_hb_into_pair_energies( True )
        self.scorefxn.set_energy_method_options( emopts )

        # randomize our root residue option
        self.randomize_root = randomize_root

        # Time test
        self.time_test = time_test

        # Debug class
        self.DEBUG = debug

        # Empty Pose test
        self.empty_pose_test = empty_pose_test


    def residues_to_perturb(self, peptide_length: int, root: int) -> List[int]:
        """
        Determine which residues should be perturbed, based on macrocycle size

        PARAMS
        ------
        :peptide_length: Amount of amino acids in our peptide
        :root: The root residue idx

        RETURNS
        -------
        :perturb_residues: A list of residues that are not pivots
        :pivote_res: A list of pivot residues
        """
        # init our output
        perturb_residues = []
        pivot_res = []
        # loop through our residues and get rid of +1, -1 from root
        for ir in range(1, peptide_length+1):
            if (ir != root+1) and (ir != root-1) and (ir != root) and (ir != peptide_length):
                perturb_residues.append(ir)
            elif ir != root:
                pivot_res.append(ir)
        # Fix Pivot RES order
        while pivot_res[0] < root:
            pre = pivot_res.pop(0)
            pivot_res.append(pre)
        return perturb_residues, pivot_res

    def get_nonroot_residues(self, peptide_length: int, root: int) -> List[int]:
        """
        Grab the residue indexes in the order that they need to be in for GenKIC as in
        from root -> peptide_length -> 1 -> root (where root is excluded)

        PARAMS
        ------
        :peptide_length: Amount of amino acids in our peptide
        :root: The root residue idx

        RETURNS
        -------
        :non_root: A list of residues that are not the root residue. Also,
            in the order that needs to be specified for genkic
        """
        # init list
        non_root = [i for i in range(root+1, peptide_length+1)] # post root
        pre_root = [i for i in range(1, root)]
        non_root.extend(pre_root)
        return non_root

    def apply_hbond_filter(self, pose: io.Pose) -> Tuple[bool, int]:
        """
        Count internal bb-bb filter where we get if it meets our cutoff

        PARAMS
        ------
        :pose: Glycine minimized pose

        RETURNS
        -------
        :filter_check: The bool output of if our peptide passes
        :hbond_count: The number of bb-bb hbonds interaction
        """
        # Get the length of our pose
        pep_len = pose.total_residue()
        root = self.foldtree_define(pep_len)

        ### init variables
        cutoff_hbonds = pep_len // 3  # we want this many internal hbonds (6: 1, 7:2, 8:2, 9:3)

        ### Setup Filters for our backbones
        # Number of internal Hbonds
        hbond_filter = PeptideInternalHbondsFilter()
        hbond_filter.set_hbond_cutoff(cutoff_hbonds)
        hbond_filter.set_hbond_types(
            backbone_backbone_setting = True,
            backbone_sidechain_setting = False,
            sidechain_sidechain_setting = False,
        )
        hbond_filter.set_scorefxn(self.scorefxn)

        # bool out
        filter_check = hbond_filter.apply(pose)
        hbond_count = hbond_filter.report_sm(pose)
        pose.scores["bb-hbonds"] = hbond_count
        return filter_check

    def anchor_randomizebyrama(self,
                               pose: io.Pose,
                               anchor_resi: int,
                               ) -> protocols.backbone_moves.RandomizeBBByRamaPrePro:
        """
        Randomize Anchor/root position by Rama Prepro Mover, since GenKIC does not alter
        it by default

        PARAMS
        ------
        :pose: Our Pose with the residue to randomize
        :anchor_resi: The residue index position of our anchor (1-indexing)
        """
        randomizeBB = protocols.backbone_moves.RandomizeBBByRamaPrePro()
        randomizeBB.set_residue_selector(
            select.ResidueIndexSelector(anchor_resi)
        )
        randomizeBB.apply(pose)


    def apply_genkic(self, pose: io.Pose) -> io.Pose:
        """
        Apply the generalized kinematic loop closure to a pose. This will designate the appropriate
        size genkic to apply

        PARAMS
        ------
        :pose: Our input glycine pose that has been set up

        RETURNS
        -------
        :genkic_pose: A stochastically sampled backbone given sequence RAMA preferences
        """
        # Get the length of our pose
        pep_len = pose.total_residue()
        root = self.foldtree_define(pep_len)

        # Calculate which residues to perturb and set as pivots
        free_residues, pivot_res = self.residues_to_perturb(pep_len, root)

        # Calculate residues to include in GenKIC
        non_root_residues = self.get_nonroot_residues(pep_len, root)
        # init the genkic class object
        GenKIC = genkic.GeneralizedKIC()
        GenKIC.set_closure_attempts(500)
        GenKIC.set_min_solution_count(1)
        GenKIC.set_selector_type("lowest_energy_selector")
        GenKIC.set_selector_scorefunction(self.scorefxn)
        # Add bb randomization for Anchor (rama prepro) if doing selection
        if self.randomize_root:
            if self.DEBUG: print("RANDOMIZE ROOT RESIDUE (THIS IS ONLY DONE FOR IN SOLUTION GENERATION)\nApplying To Pose Brefore GenKIC Closure...")
            self.anchor_randomizebyrama(pose, root)

        for ir in non_root_residues:
            GenKIC.add_loop_residue(ir)
        GenKIC.close_bond(
            rsd1 = pose.total_residue(),
            at1 = "C",
            rsd2 = 1,
            at2 = "N",
            bondlength = 1.32,
            bondangle1 = 114.,
            bondangle2 = 123.,
            torsion = 180.,
            rsd1_before = 0,
            at1_before = "",
            rsd2_after = 0,
            at2_after = "",
            randomize_this_torsion = False,
            randomize_flanking_torsions = False,
        )
        GenKIC.set_pivot_atoms(
            rsd1 = pivot_res[0],
            at1 = "CA",
            rsd2 = pivot_res[1],
            at2 = "CA",
            rsd3 = pivot_res[2],
            at3 = "CA",
        )
        GenKIC.add_perturber(genkic.perturber.perturber_effect.randomize_backbone_by_rama_prepro)
        for ir in free_residues:
            GenKIC.add_residue_to_perturber_residue_list(ir)
        GenKIC.add_filter(genkic.filter.filter_type.loop_bump_check)
        for ir in non_root_residues:
            GenKIC.add_filter(genkic.filter.filter_type.rama_prepro_check)
            GenKIC.set_filter_resnum(ir)
            GenKIC.set_filter_rama_cutoff_energy(2)

        GenKIC.apply(pose)

        return pose.clone()

    def declare_terminal_bond(self, pose: io.Pose) -> int:
       """
       Fix terminal bond

       PARAMS
       ------
       :pose: Pose object
       """
       # Fix termini, though this isn't that important
       declarebond = sm.DeclareBond()
       declarebond.set(
           res1 = pose.total_residue(),
           atom1 = 'C',
           res2 = 1,
           atom2 ='N',
           add_termini = True,
       )
       declarebond.apply(pose)
       pose.update_residue_neighbors()
       return 0

    def foldtree_define(self, N: int) -> int:
        """
        Define what the root of our fold tree is

        PARAMS
        ------
        :N (int): Desired size of input structure

        RETURNS
        -------
        :root (int): The root residue index (zero indexed)
        """
        if N % 2 == 0:
                return int(N/2)
        else:
            return int((N-1)/2)

    def full_cyclize_pose_withoutbond(self, pose: io.Pose) -> int:
        """
        Cyclize the pose using a Rosetta based method

        PARAMS
        ------
        :pose: Un cyclized pose
        """
        # Modify the variant types
        modifyvariant_nterm = sm.ModifyVariantTypeMover()
        modifyvariant_nterm.set_additional_type_to_add("CUTPOINT_UPPER")
        modifyvariant_nterm.set_residue_selector(select.ResidueIndexSelector(1))
        modifyvariant_nterm.apply(pose)

        modifyvariant_cterm = sm.ModifyVariantTypeMover()
        modifyvariant_cterm.set_additional_type_to_add("CUTPOINT_LOWER")
        modifyvariant_cterm.set_residue_selector(select.ResidueIndexSelector(pose.total_residue()))
        modifyvariant_cterm.apply(pose)

        pose.update_residue_neighbors()
        return 0

    def full_cyclize_native(self, pose: io.Pose) -> int:
        """
        Cyclize the native PDB pose using a Rosetta based method

        PARAMS
        ------
        :pose: Un cyclized pose
        """
        # Modify the variant types
        modifyvariant_nterm = sm.ModifyVariantTypeMover()
        modifyvariant_nterm.set_additional_type_to_remove("CUTPOINT_UPPER")
        modifyvariant_nterm.set_additional_type_to_remove("LOWER_TERMINUS_VARIANT")
        modifyvariant_nterm.set_additional_type_to_remove("ACETYLATED_NTERMINUS_VARIANT")
        modifyvariant_nterm.set_additional_type_to_remove("ACETYLATED_NTERMINUS_CONNECTION_VARIANT")
        modifyvariant_nterm.set_residue_selector(select.ResidueIndexSelector(1))
        modifyvariant_nterm.apply(pose)

        modifyvariant_cterm = sm.ModifyVariantTypeMover()
        modifyvariant_cterm.set_additional_type_to_remove("CUTPOINT_LOWER")
        modifyvariant_cterm.set_additional_type_to_remove("UPPER_TERMINUS_VARIANT")
        modifyvariant_cterm.set_additional_type_to_remove("DIMETHYLATED_CTERMINUS_VARIANT")
        modifyvariant_cterm.set_residue_selector(select.ResidueIndexSelector(pose.total_residue()))
        modifyvariant_cterm.apply(pose)

        # Fix termini, though this isn't that important
        declarebond = sm.DeclareBond()
        declarebond.set(
            res1 = pose.total_residue(),
            atom1 = 'C',
            res2 = 1,
            atom2 ='N',
            add_termini = True,
        )
        declarebond.apply(pose)
        return 0

    def generate_initial_empty_alanine_pose(self, N: int) -> core.pose.Pose:
       """
       Generate a single cyclic peptide pose of alanines that is of size N

       PARAMS
       ------
       :N (int): Desired size of input structure

       RETURNS
       -------
       :gly_pose (core.pose.Pose): A Rosetta Pose instance that is all alanines of desired N-mer size
           FoldTree has been updated and is uses the middle residue as a root
           Cutpoints are added to the ends and a peptide bond is declared.
       """

       if self.empty_pose_test:
           empty_pose = core.pose.Pose()
           ala_pose = peptide_stub_mover(
                   pose = empty_pose,
                   residue_type = "ALA",
                   sequence_len = N,
                   )
           flip_trans_omegabonds(ala_pose, N)
       else:
           ala_pose = io.pose_from_sequence(
               seq="A"*N,
               res_type="fa_standard",
               auto_termini=False,
           )

       # Set the fold tree
       root = self.foldtree_define(N)
       ft = core.kinematics.FoldTree()
       ft.clear()
       ft.add_edge(root, 1 -1)
       ft.add_edge(root, N, -1)
       ft.reorder(root)
       ala_pose.fold_tree(ft)

       # Cyclize the pose (but don't declare a bond yet)
       self.full_cyclize_pose_withoutbond(ala_pose)
       return ala_pose
    def generate_initial_empty_glycine_pose(self, N: int) -> core.pose.Pose:
       """
       Generate a single cyclic peptide pose of glycines that is of size N

       PARAMS
       ------
       :N (int): Desired size of input structure

       RETURNS
       -------
       :gly_pose (core.pose.Pose): A Rosetta Pose instance that is all glycines of desired N-mer size
           FoldTree has been updated and is uses the middle residue as a root
           Cutpoints are added to the ends and a peptide bond is declared.
       """

       if self.empty_pose_test:
           empty_pose = core.pose.Pose()
           gly_pose = peptide_stub_mover(
                   pose = empty_pose,
                   residue_type = "GLY",
                   sequence_len = N,
                   )
           flip_trans_omegabonds(gly_pose, N)
       else:
           gly_pose = io.pose_from_sequence(
               seq="G"*N,
               res_type="fa_standard",
               auto_termini=False,
           )

       # Set the fold tree
       root = self.foldtree_define(N)
       ft = core.kinematics.FoldTree()
       ft.clear()
       ft.add_edge(root, 1 -1)
       ft.add_edge(root, N, -1)
       ft.reorder(root)
       gly_pose.fold_tree(ft)

       # Cyclize the pose (but don't declare a bond yet)
       self.full_cyclize_pose_withoutbond(gly_pose)
       return gly_pose

    def minimize_pose(
        self,
        pose: io.Pose,
        ) -> int:
        """
        Apply a minmover to the pose for scoring

        PARAMS
        ------
        :pose: Pose with desired motif residues and the scaffold is poly-glycine

        RETURNS
        -------
        :pose: minimized version of this pose
        """
        pose_before = pose.clone()

        # Set up a movemap
        mm = core.kinematics.MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)

        # Set up minimover
        minmover = protocols.minimization_packing.MinMover()
        minmover.movemap(mm)
        minmover.score_function(self.scorefxn)
        minmover.min_type("lbfgs_armijo_nonmonotone")
        minmover.tolerance(1.0e-7)
        minmover.apply(pose)

        rmsd = core.scoring.bb_rmsd_including_O(pose_before, pose)
        pose.scores["MinMover_bbheavy_RMSD"] = rmsd
        self.scorefxn(pose)
        return 0

    @timeit
    def generate_ensemble(self, s: int, nstruct: int, nofilter: bool, alanine: bool) -> int:
        """
        Helper function for running full process of pose gen + genkic/filter + minimize output
        for a given size (s)

        PARAMS
        ------
        :s: The desired output size of our pose
        :nstruct: Number of desired outputs
        :nofilter: Argparse argument for if you should filter or not based on hbonds
        :alanine: True if our pose are going to be made of alanines, this is to make
            all L-alpha-amino acid conformations

        RETURNS
        -------
        :ensemble: No actual return, but instead a silent file of starting poses for design
        """
        print("-"*4, "GenKIC, Size:", s, "Nstruct:", nstruct, "DEBUG:", self.DEBUG, "Alanine:", alanine, "-"*4)

        # Generate our initial pose
        if alanine:
            pose = self.generate_initial_empty_alanine_pose(s)
        else:
            pose = self.generate_initial_empty_glycine_pose(s)

        # Declare our terminal bond
        self.declare_terminal_bond(pose)

        # Setup SilentFile Output
        opts = silent.SilentFileOptions()
        opts.in_fullatom(True)
        opts.set_binary_output(True)
        silentFile = silent.SilentFileData(opts)
        if self.time_test:
            out_name = f"genkicbb_size{s}_{datetime.date.today().strftime('%m%d%Y')}_timetest.silent"
        else:
            out_name = f"genkicbb_size{s}_{datetime.date.today().strftime('%m%d%Y')}.silent"

        # Apply genkic to our pose
        success = 0
        while success < nstruct:
        # for n in range(nstruct):
            # We clone our pose here, so that the bb torsion randomization is stochastic
            pre_genkic_pose = pose.clone()
            # Apply Genkic to pose and clone it
            genkic_pose = self.apply_genkic(pre_genkic_pose)
            # Minimize with MinMover for small energetic improvement
            self.minimize_pose(genkic_pose)
            # Calculate bb-bb hbond amount + Get True/False output
            hbond_filter = self.apply_hbond_filter(genkic_pose)
            # We want to filter by bb-bb hbonds (HbondFilter: True + use-hbond-filter: True)
            if hbond_filter and not nofilter:
                # Write to silent file
                silentStruct = silent.BinarySilentStruct(
                    opts, genkic_pose, f"size{s}_bb_{str(success+1).zfill(6)}",
                )
                if self.DEBUG:
                    genkic_pose.dump_pdb(f"debug_{s}_{str(success+1).zfill(6)}.pdb")
                success+=1
                if success % 1000 == 0:
                    print("Successfully Generated Filtered Poly-glycine Backbones:", success)
                if self.time_test:
                    silentFile.write_silent_struct(silentStruct, out_name)
                else:
                    silentFile.add_structure(silentStruct)
                # We don't want to filter by number of bb-bb hbonds. So Flag not set use-hbond-filter: False
            elif nofilter:
                # Write to silent file
                silentStruct = silent.BinarySilentStruct(
                    opts, genkic_pose, f"size{s}_bb_{str(success+1).zfill(6)}",
                )
                if self.DEBUG:
                    genkic_pose.dump_pdb(f"debug_{s}_{str(success+1).zfill(6)}.pdb")
                success+=1
                if success % 1000 == 0:
                    print("Successfully Generated Non-Filtered Poly-glycine Backbones:", success)
                if self.time_test:
                    silentFile.write_silent_struct(silentStruct, out_name)
                else:
                    silentFile.add_structure(silentStruct)

        if not self.time_test:
            silentFile.write_all(out_name)
        print("SilentFile Finished:", out_name)

        return 0

if __name__ == "__main__":
    p = argparse.ArgumentParser("Generate Macrocycle Backbones",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                )
    p.add_argument('-s', '--size', nargs="+", type=int, default=[6, 7, 8, 9], help="List of number of \
            residue macrocycles you'd like to generate. (e.g. --size 6 7)")
    p.add_argument('-n', '--nstruct', type=int, default=10_000, help="Number of structures to generate per size.")
    p.add_argument('--debug', action="store_true", help="Dump PDBs for testing, unmute Rosetta, print helpful trace messages")
    p.add_argument('--nofilter', action="store_true", help="Dont apply strict filter on hbonds")
    p.add_argument('--sample-root', action="store_true", help="If generating samples in-solution then this should be set \
            as your root residue is generally not an anchor, which means it can sample freely..")
    p.add_argument("--time-test", action="store_true", help="This is only for a time comparison between XML and PyRosetta \
            as XML writes to disk for every structure, but here we do not. This can be done for time comparisons. Dont set \
            as it will make the process slower.")
    p.add_argument("--empty-pose-test", action="store_true", help="Run an empty pose test, since this is more representative of the XML script")
    p.add_argument("--alanine", action="store_true", help="Make our pose all Alanines, this is important if you want your conformations to be \
            all L-alpha-amino acids, if not passed then pose will be all Glycines. Which allow heterochiral conformations."
    args = p.parse_args()

    bbgen = BackboneGeneration(args.debug, args.sample_root, args.time_test, args.empty_pose_test)
    for s in args.size:
        bbgen.generate_ensemble(s, args.nstruct, args.nofilter, args.alanine)

