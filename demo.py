


# |%%--%%| <FHk2C2mAci|tveJZTmZdq>
r"""°°°
# QCArchive+QCMLForge Demo with CyberShuttle

The first half of this demo shows how to use QCArchive to setup a dataset
and run computations with ease. The compute resource can be local or
through CyberShuttle. 

The second half of this demo shows how one can consume the generated data
to train AP-Net2 and dAPNet2 models through QCMLForge. 
°°°"""
# |%%--%%| <tveJZTmZdq|sb2BSlStsm>

import psi4
from pprint import pprint as pp
import pandas as pd
import numpy as np
import re
from qm_tools_aw import tools
from pprint import pprint as pp
# QCElemental Imports
from qcelemental.models import Molecule
import qcelemental as qcel
# Dataset Imports
from qcportal import PortalClient
from qcportal.singlepoint import SinglepointDataset, SinglepointDatasetEntry, QCSpecification
from qcportal.manybody import ManybodyDataset, ManybodyDatasetEntry, ManybodyDatasetSpecification, ManybodySpecification
from torch import manual_seed

manual_seed(42)

h2kcalmol = qcel.constants.hartree2kcalmol

# |%%--%%| <sb2BSlStsm|78H3oHPXBB>
r"""°°°
# QCArchive Setup
°°°"""
# |%%--%%| <78H3oHPXBB|BVc6W6uOta>

from setup_qcfractal import setup_qcarchive_qcfractal
import os

setup_qcarchive_qcfractal(
    QCF_BASE_FOLDER=os.path.join(os.getcwd(), "qcfractal"),
    start=False,
    reset=False,
    db_config={
        "name": None,
        "enable_security": "false",
        "allow_unauthenticated_read": None,
        "logfile": None,
        "loglevel": None,
        "service_frequency": 5,
        "max_active_services": None,
        "heartbeat_frequency": 60,
        "log_access": None,
        "database": {
            "base_folder": None,
            "host": None,
            "port": 5432,
            "database_name": "qca",
            "username": None,
            "password": None,
            "own": None,
        },
        "api": {
            "host": None,
            "port": 7777,
            "secret_key": None,
            "jwt_secret_key": None,
        },
    },
    resources_config={
            "update_frequency": 5,
            "cores_per_worker": 8,
            "max_workers": 3,
            "memory_per_worker": 20,
    }
)

# |%%--%%| <BVc6W6uOta|CwEhpqwLXX>

!qcfractal-server --config=`pwd`/qcfractal/qcfractal_config.yaml start > qcfractal_server.log & disown

# |%%--%%| <CwEhpqwLXX|3HjtiyIuFg>

!qcfractal-compute-manager --config=`pwd`/qcfractal/resources.yml > qcfractal_compute.log & disown

# |%%--%%| <3HjtiyIuFg|2gmdKuBAuk>

# NOTE kill server when finished by removing the # and executing:
# !ps aux | grep qcfractal | awk '{ print $2 }' | xargs kill -9

# |%%--%%| <2gmdKuBAuk|X31FadtbG5>
r"""°°°
# QCArchive single point example
°°°"""
# |%%--%%| <X31FadtbG5|hMCmRgdGJ4>

# Running a single job
client = PortalClient("http://localhost:7777", verify=False)
mol = Molecule.from_data(
    """
     0 1
     O  -1.551007  -0.114520   0.000000
     H  -1.934259   0.762503   0.000000
     H  -0.599677   0.040712   0.000000
     --
     0 1
     O   1.350625   0.111469   0.000000
     H   1.680398  -0.373741  -0.758561
     H   1.680398  -0.373741   0.758561

     units angstrom
     no_reorient
     symmetry c1
"""
)

psi4.set_options(
    {"basis": "aug-cc-pvdz", "scf_type": "df", "e_convergence": 6, "freeze_core": True}
)

client.add_singlepoints(
    [mol],
    "psi4",
    driver="energy",
    method="b3lyp",
    basis="aug-cc-pvdz",
    keywords={"scf_type": "df", "e_convergence": 6, "freeze_core": True},
    tag="local",
)

# Can print records
# for rec in client.query_records():
#     pp(rec.dict)
#     pp(rec.error)

# |%%--%%| <hMCmRgdGJ4|6DmD5wp1nB>
r"""°°°
# QCArchive dataset examples
°°°"""
# |%%--%%| <6DmD5wp1nB|Z0wXrcgRq8>

# Creating a QCArchive Dataset...
# Load in a dataset from a recent Sherrill work (Levels of SAPT II)
df_LoS = pd.read_pickle("./combined_df_subset_358.pkl")
print(df_LoS[['Benchmark', 'SAPT2+3(CCD)DMP2 TOTAL ENERGY aqz', 'MP2 IE atz', 'SAPT0 TOTAL ENERGY adz' ]])

# Limit to 100 molecules with maximum of 16 atoms to keep computational cost down
df_LoS['size'] = df_LoS['atomic_numbers'].apply(lambda x: len(x))
df_LoS = df_LoS[df_LoS['size'] <= 16]
df_LoS = df_LoS.sample(200, random_state=42, axis=0).copy()
df_LoS.reset_index(drop=True, inplace=True)
print(df_LoS['size'].describe())

# Create QCElemntal Molecules to generate the dataset
def qcel_mols(row):
    """
    Convert the row to a qcel molecule
    """
    atomic_numbers = [row['atomic_numbers'][row['monAs']], row['atomic_numbers'][row['monBs']]]
    coords = [row['coordinates'][row['monAs']], row['coordinates'][row['monBs']]]
    cm = [
        [row['monA_charge'], row['monA_multiplicity']],
        [row['monB_charge'], row['monB_multiplicity']],
     ]
    return tools.convert_pos_carts_to_mol(atomic_numbers, coords, cm)
df_LoS['qcel_molecule'] = df_LoS.apply(qcel_mols, axis=1)
geoms = df_LoS['qcel_molecule'].tolist()
ref_IEs = df_LoS['Benchmark'].tolist()
sapt0_adz = (df_LoS['SAPT0 TOTAL ENERGY adz'] * h2kcalmol).tolist()

# |%%--%%| <Z0wXrcgRq8|iwmcvViziS>
r"""°°°
## Singlepoint Dataset
°°°"""
# |%%--%%| <iwmcvViziS|i8ICwzPWaD>

# Create client dataset

ds_name = 'S22-singlepoint'
client_datasets = [i['dataset_name'] for i in client.list_datasets()]
# Check if dataset already exists, if not create a new one
if ds_name not in client_datasets:
    ds = client.add_dataset("singlepoint", ds_name,
                            f"Dataset to contain {ds_name}")
    print(f"Added {ds_name} as dataset")
    # Insert entries into dataset
    entry_list = []
    for idx, mol in enumerate(geoms):
        extras = {
            "name": 'S22-' + str(idx),
            "idx": idx,
        }
        mol = Molecule.from_data(mol.dict(), extras=extras)
        ent = SinglepointDatasetEntry(name=extras['name'], molecule=mol)
        entry_list.append(ent)
    ds.add_entries(entry_list)
    print(f"Added {len(entry_list)} molecules to dataset")
else:
    ds = client.get_dataset("singlepoint", ds_name)
    print(f"Found {ds_name} dataset, using this instead")

print(ds)

# |%%--%%| <i8ICwzPWaD|7ZAOPlzuUX>

# Can delete the dataset if you want to start over. Need to know dataset_id
# client.delete_dataset(dataset_id=ds.id, delete_records=True)

# |%%--%%| <7ZAOPlzuUX|dSW1A9HxYB>

# Multipole Example
# method, basis = "hf", "sto-3g"
#
# # Set the QCSpecification (QM interaction energy in our case)
# spec = QCSpecification(
#     program="psi4",
#     driver="energy",
#     method=method,
#     basis=basis,
#     keywords={
#         "d_convergence": 8,
#         "dft_radial_points": 99,
#         "dft_spherical_points": 590,
#         "e_convergence": 10,
#         "guess": "sad",
#         "mbis_d_convergence": 9,
#         "mbis_radial_points": 99,
#         "mbis_spherical_points": 590,
#         "scf_properties": ["mbis_charges", "MBIS_VOLUME_RATIOS"],
#         "scf_type": "df",
#     },
#     protocols={"wavefunction": "orbitals_and_eigenvalues"},
# )
# ds.add_specification(name=f"psi4/{method}/{basis}", specification=spec)

# |%%--%%| <dSW1A9HxYB|vMkm00fo00>

# SAPT0 Example
method, basis = "SAPT0", "cc-pvdz"

# Set the QCSpecification (QM interaction energy in our case)
spec = QCSpecification(
    program="psi4",
    driver="energy",
    method=method,
    basis=basis,
    keywords={
        "scf_type": "df",
    },
)
ds.add_specification(name=f"psi4/{method}/{basis}", specification=spec)

# |%%--%%| <vMkm00fo00|dwYb9dbQNI>

# Run the computations
ds.submit()
print(f"Submitted {ds_name} dataset")

# |%%--%%| <dwYb9dbQNI|2JMCNlehez>

# Check the status of the dataset - can repeatedly run this to see the progress
ds.status()

# |%%--%%| <2JMCNlehez|kT14mJyNqJ>
r"""°°°
## Manybody Dataset
°°°"""
# |%%--%%| <kT14mJyNqJ|g31JlHrgso>

# Create client dataset
ds_name_mb = 'S22-manybody'
client_datasets = [i['dataset_name'] for i in client.list_datasets()]
# Check if dataset already exists, if not create a new one
if ds_name_mb not in client_datasets:
    print("Setting up new dataset:", ds_name_mb)
    ds_mb = client.add_dataset("manybody", ds_name_mb,
                            f"Dataset to contain {ds_name_mb}")
    print(f"Added {ds_name_mb} as dataset")
    # Insert entries into dataset
    entry_list = []
    for idx, mol in enumerate(geoms):
        ent = ManybodyDatasetEntry(name=f"S22-IE-{idx}", initial_molecule=mol)
        entry_list.append(ent)
    ds_mb.add_entries(entry_list)
    print(f"Added {len(entry_list)} molecules to dataset")
else:
    ds_mb = client.get_dataset("manybody", ds_name_mb)
    print(f"Found {ds_name_mb} dataset, using this instead")

print(ds_mb)

# Can delete the dataset if you want to start over. Need to know dataset_id
# client.delete_dataset(dataset_id=2, delete_records=True)

# |%%--%%| <g31JlHrgso|bYERcUudd0>

ds_mb.status()

# |%%--%%| <bYERcUudd0|gauw3VIjl9>

# Set multiple levels of theory - you can add/remove levels as you desire.
# Computational scaling will get quite expensive with better methods and larger
# basis sets

methods = [
    'hf', 'pbe',
]
basis_sets = [
    '6-31g*'
]

for method in methods:
    for basis in basis_sets:
        # Set the QCSpecification (QM interaction energy in our case)
        qc_spec_mb = QCSpecification(
            program="psi4",
            driver="energy",
            method=method,
            basis=basis,
            keywords={
                "d_convergence": 8,
                "scf_type": "df",
            },
        )

        spec_mb = ManybodySpecification(
            program='qcmanybody',
            bsse_correction=['cp'],
            levels={
                1: qc_spec_mb,
                2: qc_spec_mb,
            },
        )
        print("spec_mb", spec_mb)

        ds_mb.add_specification(name=f"psi4/{method}/{basis}", specification=spec_mb)

        # Run the computations
        ds_mb.submit()
        print(f"Submitted {ds_name} dataset")
# Check the status of the dataset - can repeatedly run this to see the progress
ds_mb.status()

# |%%--%%| <gauw3VIjl9|qYukdPBXmi>

pp(ds.status())
pp(ds_mb.status())

# |%%--%%| <qYukdPBXmi|ChCOdcBiXj>

pp(ds)
pp(ds_mb)
pp(ds_mb.computed_properties)

# |%%--%%| <ChCOdcBiXj|xgdzc0Klhx>
r"""°°°
# Data Assembly

While you can execute the following blocks before all computations are complete, it is recommended to wait until all computations are complete to continue.
°°°"""
# |%%--%%| <xgdzc0Klhx|7RHL31QOoC>

# Singlepoint data assemble
def assemble_singlepoint_data(record):
    record_dict = record.dict()
    qcvars = record_dict["properties"]
    level_of_theory = f"{record_dict['specification']['method']}/{record_dict['specification']['basis']}"
    sapt_energies = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
    if "mbis charges" in qcvars:
        charges = qcvars["mbis charges"]
        dipoles = qcvars["mbis dipoles"]
        quadrupoles = qcvars["mbis quadrupoles"]
        n = len(charges)
        charges = np.reshape(charges, (n, 1))
        dipoles = np.reshape(dipoles, (n, 3))
        quad = np.reshape(quadrupoles, (n, 3, 3))

        quad = [q[np.triu_indices(3)] for q in quad]
        quadrupoles = np.array(quad)
        multipoles = np.concatenate(
            [charges, dipoles, quadrupoles], axis=1)
        return (
        record.molecule,
        qcvars['mbis volume ratios'],
        qcvars['mbis valence widths'],
        qcvars['mbis radial moments <r^2>'],
        qcvars['mbis radial moments <r^3>'],
        qcvars['mbis radial moments <r^4>'],
        record.molecule.atomic_numbers,
        record.molecule.geometry * qcel.constants.bohr2angstroms,
        multipoles,
        int(record.molecule.molecular_charge),
        record.molecule.molecular_multiplicity,
        sapt_energies,
        )
    else:
        sapt_energies[0] = qcvars['sapt total energy']
        sapt_energies[1] = qcvars['sapt elst energy']
        sapt_energies[2] = qcvars['sapt exch energy']
        sapt_energies[3] = qcvars['sapt ind energy']
        sapt_energies[4] = qcvars['sapt disp energy']
        return (
        record.molecule,
        None,
        None,
        None,
        None,
        None,
        record.molecule.atomic_numbers,
        record.molecule.geometry * qcel.constants.bohr2angstroms,
        None,
        int(record.molecule.molecular_charge),
        record.molecule.molecular_multiplicity,
        sapt_energies,
        )




def assemble_singlepoint_data_value_names():
    return [
        'qcel_molecule',
        "volume ratios",
        "valence widths",
        "radial moments <r^2>",
        "radial moments <r^3>",
        "radial moments <r^4>",
        "Z",
        "R",
        "cartesian_multipoles",
        "TQ",
        "molecular_multiplicity",
        "SAPT Energies"
    ]

df = ds.compile_values(
    value_call=assemble_singlepoint_data,
    value_names=assemble_singlepoint_data_value_names(),
    unpack=True,
)
pp(df.columns.tolist())
df_sapt0 = df['psi4/SAPT0/cc-pvdz']

# |%%--%%| <7RHL31QOoC|CDD1QjxHpc>

def assemble_data(record):
    record_dict = record.dict()
    qcvars = record_dict["properties"]
    level_of_theory = f"{record_dict['specification']['levels'][2]['method']}/{record_dict['specification']['levels'][2]['basis']}"
    CP_IE = qcvars['results']['cp_corrected_interaction_energy'] * h2kcalmol
    NOCP_IE = qcvars['results'].get('nocp_corrected_interaction_energy', np.nan) * h2kcalmol
    return (
    record.initial_molecule,
    CP_IE,
    NOCP_IE,
    record.initial_molecule.atomic_numbers,
    record.initial_molecule.geometry * qcel.constants.bohr2angstroms,
    int(record.initial_molecule.molecular_charge),
    record.initial_molecule.molecular_multiplicity,
    )

def assemble_data_value_names():
    return [
        'qcel_molecule',
        "CP_IE",
        "NOCP_IE",
        "Z",
        "R",
        "TQ",
        "molecular_multiplicity"
    ]

df_mb = ds_mb.compile_values(
    value_call=assemble_data,
    value_names=assemble_data_value_names(),
    unpack=True,
)

pp(df_mb.columns.tolist())

# |%%--%%| <CDD1QjxHpc|XT87RegBfm>

from cdsg_plot import error_statistics

df_sapt0['sapt0 total energes'] = df_sapt0['SAPT Energies'].apply(lambda x: x[0] * h2kcalmol)
df_plot = pd.DataFrame(
    {
        "qcel_molecule": df_mb["psi4/pbe/6-31g*"]["qcel_molecule"],
        "HF/6-31G*": df_mb["psi4/hf/6-31g*"]["CP_IE"],
        "PBE/6-31G*": df_mb["psi4/pbe/6-31g*"]["CP_IE"],
        'SAPT0/cc-pvdz': df_sapt0['sapt0 total energes'].values,
    }
)
print(df_plot)
id = [int(i[7:]) for i in df_plot.index]
df_plot['id'] = id
df_plot.sort_values(by='id', inplace=True, ascending=True)
df_plot['reference'] = ref_IEs
df_plot['SAPT0/aug-cc-pvdz'] = sapt0_adz
df_plot['HF/6-31G* error'] = (df_plot['HF/6-31G*'] - df_plot['reference']).astype(float)
df_plot['PBE/6-31G* error'] = (df_plot['PBE/6-31G*'] - df_plot['reference']).astype(float)
df_plot['SAPT0/cc-pvdz error'] = (df_plot['SAPT0/cc-pvdz'] - df_plot['reference']).astype(float)
df_plot['SAPT0/aug-cc-pvdz error'] = (df_plot['SAPT0/aug-cc-pvdz'] - df_plot['reference']).astype(float)
pd.set_option('display.max_rows', None)
print(df_plot[['PBE/6-31G*', 'SAPT0/cc-pvdz', 'reference', "SAPT0/aug-cc-pvdz"]])
print(df_plot[['PBE/6-31G* error', 'SAPT0/cc-pvdz error', "SAPT0/aug-cc-pvdz error"]].describe())

# |%%--%%| <XT87RegBfm|6ziYQ8PdtF>
r"""°°°
# Plotting the interaction energy errors
°°°"""
# |%%--%%| <6ziYQ8PdtF|dRiuyCtOh1>

error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        # "B3LYP/6-31G*": "B3LYP/6-31G* error",
        "SAPT0/cc-pvdz": "SAPT0/cc-pvdz error",
        "SAPT0/aug-cc-pvdz": "SAPT0/aug-cc-pvdz error",
    },
    output_filename="S22-IE.png",
    figure_size=(6, 6),
    x_label_fontsize=16,
    ylim=(-15, 15),
    rcParams={},
    usetex=False,
    ylabel=r"IE Error vs. CCSD(T)/CBS (kcal/mol)",
)

# |%%--%%| <dRiuyCtOh1|NoxyyvkpUK>
r"""°°°
# QCMLForge

## AP-Net2 inference
°°°"""
# |%%--%%| <NoxyyvkpUK|8jtfD3m3S4>

import apnet_pt
from apnet_pt.AtomPairwiseModels.apnet2 import APNet2Model
from apnet_pt.AtomModels.ap2_atom_model import AtomModel

atom_model = AtomModel().set_pretrained_model(model_id=0)
ap2 = APNet2Model(atom_model=atom_model.model).set_pretrained_model(model_id=0)
ap2.atom_model = atom_model.model
apnet2_ies_predicted = ap2.predict_qcel_mols(
    mols=df_plot['qcel_molecule'].tolist(),
    batch_size=16
)

# |%%--%%| <8jtfD3m3S4|nxmLPrfAfx>

# AP-Net2 IE
df_plot['APNet2'] = np.sum(apnet2_ies_predicted, axis=1)
df_plot['APNet2 error'] = (df_plot['APNet2'] - df_plot['reference']).astype(float)
print(df_plot.sort_values(by='APNet2 error', ascending=True)[['APNet2', 'reference']])
error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "SAPT0/cc-pvdz": "SAPT0/cc-pvdz error",
        "SAPT0/aug-cc-pvdz": "SAPT0/aug-cc-pvdz error",
        "APNet2": "APNet2 error",
    },
    output_filename="S22-IE-AP2.png",
    rcParams={},
    usetex=False,
    figure_size=(4, 4),
    ylabel=r"IE Error vs. CCSD(T)/CBS (kcal/mol)",
)

# |%%--%%| <nxmLPrfAfx|XtMJQEokjd>

# Training models on new QM data: Transfer Learning

from apnet_pt import pairwise_datasets

ds2 = pairwise_datasets.apnet2_module_dataset(
    root="data_dir",
    spec_type=None,
    atom_model=atom_model,
    qcel_molecules=df_plot['qcel_molecule'].tolist(),
    energy_labels=[np.array([i]) for i in df_plot['reference'].tolist()],
    skip_compile=True,
    force_reprocess=True,
    atomic_batch_size=8,
    prebatched=False,
    in_memory=True,
    batch_size=4,
)
print(ds2)

# |%%--%%| <XtMJQEokjd|4wjE6QG52G>
r"""°°°
## Transfer Learning
°°°"""
# |%%--%%| <4wjE6QG52G|xXslqNQSRI>

# Transfer Learning APNet2 model on computed QM data
ap2.train(
    dataset=ds2,
    n_epochs=50,
    transfer_learning=True,
    skip_compile=True,
    model_path="apnet2_transfer_learning.pt",
    split_percent=0.8,
)

# |%%--%%| <xXslqNQSRI|KPOeFMBhWm>

# AP-Net2 IE
apnet2_ies_predicted_transfer = ap2.predict_qcel_mols(
    mols=df_plot['qcel_molecule'].tolist(),
    batch_size=16,
)
df_plot['APNet2 transfer'] = np.sum(apnet2_ies_predicted_transfer, axis=1)
df_plot['APNet2 transfer error'] = (df_plot['APNet2 transfer'] - df_plot['reference']).astype(float)

error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "SAPT0/aug-cc-pvdz": "SAPT0/aug-cc-pvdz error",
        "APNet2": "APNet2 error",
        "APNet2 transfer": "APNet2 transfer error",
    },
    output_filename="S22-IE-AP2.png",
    rcParams={},
    usetex=False,
    figure_size=(6, 4),
    ylabel=r"IE Error vs. CCSD(T)/CBS (kcal/mol)",
)

# |%%--%%| <KPOeFMBhWm|M3Q6tUC2AJ>
r"""°°°
## $\Delta$AP-Net2
°°°"""
# |%%--%%| <M3Q6tUC2AJ|ROLVxSHnj2>

from apnet_pt.pt_datasets.dapnet_ds import dapnet2_module_dataset_apnetStored

delta_energies = df_plot['HF/6-31G* error'].tolist()

# Only operates in pre-batched mode
ds = dapnet2_module_dataset_apnetStored(
    root="data_dir",
    r_cut=5.0,
    r_cut_im=8.0,
    spec_type=None,
    max_size=None,
    force_reprocess=True,
    batch_size=2,
    num_devices=1,
    skip_processed=False,
    skip_compile=True,
    print_level=2,
    in_memory=True,
    m1="HF/6-31G*",
    m2="CCSD(T)/CBS",
    qcel_molecules=df_plot['qcel_molecule'].tolist(),
    energy_labels=delta_energies,
)
print(ds)

# |%%--%%| <ROLVxSHnj2|jfEVLnbL5w>

from apnet_pt.AtomPairwiseModels.dapnet2 import dAPNet2Model

dap2 = dAPNet2Model(
    atom_model=AtomModel().set_pretrained_model(model_id=0),
    apnet2_model=APNet2Model().set_pretrained_model(model_id=0).set_return_hidden_states(True),
)
dap2.train(
    ds,
    n_epochs=50,
    skip_compile=True,
    split_percent=0.6,
)

# |%%--%%| <jfEVLnbL5w|BbSVmmebsN>

dAPNet2_ies_predicted_transfer = dap2.predict_qcel_mols(
    mols=df_plot['qcel_molecule'].tolist(),
    batch_size=2,
)
df_plot['dAPNet2'] = dAPNet2_ies_predicted_transfer
df_plot['HF/6-31G*-dAPNet2'] = df_plot['HF/6-31G*'] - df_plot['dAPNet2']
print(df_plot[['dAPNet2', 'HF/6-31G*', 'HF/6-31G*-dAPNet2',  'reference']])
df_plot['dAPNet2 error'] = (df_plot['HF/6-31G*-dAPNet2'] - df_plot['reference']).astype(float)

error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "SAPT0/aug-cc-pvdz": "SAPT0/aug-cc-pvdz error",
        "APNet2": "APNet2 error",
        "APNet2 transfer": "APNet2 transfer error",
        "dAPNet2 HF/6-31G*->CCSD(T)/CBS": "dAPNet2 error",
    },
    
    rcParams={},
    usetex=False,
    figure_size=(6, 4),
    ylabel=r"IE Error vs. CCSD(T)/CBS (kcal/mol)",
)

# |%%--%%| <BbSVmmebsN|OVVcYRXWUA>

# Be careful with this for it can corrupt running status...
# !ps aux | grep qcfractal | awk '{ print $2 }' | xargs kill -9

# |%%--%%| <OVVcYRXWUA|H2TTLmA061>
r"""°°°
# The end...
°°°"""
# |%%--%%| <H2TTLmA061|ULuCE0r8ED>


