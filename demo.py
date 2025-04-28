
import psi4
from qcportal import PortalClient
from pprint import pprint as pp
from qcelemental.models import Molecule
import qcelemental as qcel
from qcportal.singlepoint import SinglepointDataset, SinglepointDatasetEntry, QCSpecification
import pandas as pd
import numpy as np
import re
from qm_tools_aw import tools

# need manybodydataset
from qcportal.manybody import ManybodyDataset, ManybodyDatasetEntry, ManybodyDatasetSpecification, ManybodySpecification
from pprint import pprint as pp

# |%%--%%| <sb2BSlStsm|BVc6W6uOta>

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

# NOTE kill server when finished by running:
#     ps aux | grep qcfractal-server | awk '{ print $2 }'
#     kill -9 <PID>

# |%%--%%| <CwEhpqwLXX|3HjtiyIuFg>

!qcfractal-compute-manager --config=`pwd`/qcfractal/resources.yml > qcfractal_compute.log & disown
# NOTE kill server when finished by running:;
#    ps aux | grep qcfractal-server | awk '{ print $2 }'
#    kill -9 <PID>

# |%%--%%| <3HjtiyIuFg|hMCmRgdGJ4>

# Running a single job
client = PortalClient("http://localhost:7777", verify=False)
# for rec in client.query_records():
#     pp(rec)

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

# for rec in client.query_records():
#     pp(rec.dict)
#     pp(rec.error)

# |%%--%%| <hMCmRgdGJ4|Z0wXrcgRq8>

# Now create S22 Interaction Energy Dataset
from s22 import geoms, ref_IEs

assert len(geoms) == len(ref_IEs), "Number of geometries and reference interaction energies do not match"

# geoms is a list of qcelemental Molecule objects that can be used to create a
# QCArchive dataset
print(len(geoms), geoms)

# |%%--%%| <Z0wXrcgRq8|i8ICwzPWaD>

# Create client dataset

ds_name = 'S22-multipoles'

try:
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
except Exception:
    ds = client.get_dataset("singlepoint", ds_name)
    print(f"Found {ds_name} dataset, using this instead")
    print(ds)

# |%%--%%| <i8ICwzPWaD|dSW1A9HxYB>

# Set the method and basis for lower requirements?
method, basis = "hf", "sto-3g"

# Set the QCSpecification (QM interaction energy in our case)
spec = QCSpecification(
    program="psi4",
    driver="energy",
    method=method,
    basis=basis,
    keywords={
        "d_convergence": 8,
        "dft_radial_points": 99,
        "dft_spherical_points": 590,
        "e_convergence": 10,
        "guess": "sad",
        "mbis_d_convergence": 9,
        "mbis_radial_points": 99,
        "mbis_spherical_points": 590,
        "scf_properties": ["mbis_charges", "MBIS_VOLUME_RATIOS"],
        "scf_type": "df",
    },
    protocols={"wavefunction": "orbitals_and_eigenvalues"},
)
ds.add_specification(name=f"psi4/{method}/{basis}", specification=spec)

# |%%--%%| <dSW1A9HxYB|vMkm00fo00>

# Set the method and basis for lower requirements?
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

# |%%--%%| <2JMCNlehez|g31JlHrgso>

# Create client dataset

ds_name_mb = 'S22-manybody'

try:
    ds_mb = client.add_dataset("manybody", ds_name_mb,
                            f"Dataset to contain {ds_name_mb}")
    print(f"Added {ds_name_mb} as dataset")
except Exception:
    ds_mb = client.get_dataset("manybody", ds_name_mb)
    print(f"Found {ds_name_mb} dataset, using this instead")
    print(ds)

# Insert entries into dataset

entry_list = []
for idx, mol in enumerate(geoms):
    print(mol)
    ent = ManybodyDatasetEntry(name=f"S22-IE-{idx}", initial_molecule=mol)
    entry_list.append(ent)
ds_mb.add_entries(entry_list)
print(f"Added {len(entry_list)} molecules to dataset")

# Set the method and basis for lower requirements?
method, basis = "hf", "sto-3g"

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
    bsse_correction=['cp', 'nocp'],
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

# |%%--%%| <g31JlHrgso|bYERcUudd0>

ds_mb.status()

# |%%--%%| <bYERcUudd0|gauw3VIjl9>

# Want multiple levels of theory

methods = [
    'hf', 'pbe', 'b3lyp',
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
# client.delete_dataset(dataset_id=2, delete_records=True)

# |%%--%%| <ChCOdcBiXj|7RHL31QOoC>

# Multipole Molecule assemble
def assemble_multipole_data(record):
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
        # pp(record_dict)
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




def assemble_multipole_data_value_names():
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
    value_call=assemble_multipole_data,
    value_names=assemble_multipole_data_value_names(),
    unpack=True,
)
# print(df)
pp(df.columns.tolist())
print(df['psi4/hf/sto-3g'])
pp(df['psi4/hf/sto-3g'].columns.tolist())
df_sapt0 = df['psi4/SAPT0/cc-pvdz']
print(df_sapt0.columns.tolist())
print(df_sapt0)

# |%%--%%| <7RHL31QOoC|CDD1QjxHpc>

h2kcalmol = qcel.constants.hartree2kcalmol

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

print(df_mb)

# |%%--%%| <CDD1QjxHpc|XT87RegBfm>

from cdsg_plot import error_statistics

pp(df_mb.columns.tolist())
print(len(df_sapt0))
print(len(df_sapt0['SAPT Energies']))
print(len(df_sapt0['SAPT Energies'].apply(lambda x: x[0] * h2kcalmol)))

df_sapt0['sapt0 total energes'] = df_sapt0['SAPT Energies'].apply(lambda x: x[0] * h2kcalmol)
print(df_sapt0['sapt0 total energes'])
df_plot = pd.DataFrame(
    {
        "qcel_molecule": df_mb["psi4/hf/sto-3g"]["qcel_molecule"],
        "HF/6-31G*": df_mb["psi4/hf/6-31g*"]["CP_IE"],
        "PBE/6-31G*": df_mb["psi4/pbe/6-31g*"]["CP_IE"],
        "B3LYP/6-31G*": df_mb["psi4/b3lyp/6-31g*"]["CP_IE"],
        'SAPT0/cc-pvdz': df_sapt0['sapt0 total energes'].values,
    }
)
print(df_plot)
id = [int(i[7:]) for i in df_plot.index]
df_plot['id'] = id
df_plot.sort_values(by='id', inplace=True, ascending=True)
df_plot['reference'] = ref_IEs
df_plot['HF/6-31G* error'] = (df_plot['HF/6-31G*'] - df_plot['reference']).astype(float)
df_plot['PBE/6-31G* error'] = (df_plot['PBE/6-31G*'] - df_plot['reference']).astype(float)
df_plot['B3LYP/6-31G* error'] = (df_plot['B3LYP/6-31G*'] - df_plot['reference']).astype(float)
df_plot['SAPT0/cc-pvdz error'] = (df_plot['SAPT0/cc-pvdz'] - df_plot['reference']).astype(float)
pd.set_option('display.max_rows', None)
print(df_plot)
print(df_plot[['HF/6-31G* error', 'PBE/6-31G* error', 'B3LYP/6-31G* error']].describe())

# |%%--%%| <XT87RegBfm|dRiuyCtOh1>

error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "B3LYP/6-31G*": "B3LYP/6-31G* error",
        "SAPT0/cc-pvdz": "SAPT0/cc-pvdz error",
    },
    output_filename="S22-IE.png",
    figure_size=(8, 6),
)

# |%%--%%| <dRiuyCtOh1|8jtfD3m3S4>

import apnet_pt
from apnet_pt.AtomPairwiseModels.apnet2 import APNet2Model

ap2 = APNet2Model().set_pretrained_model(model_id=0)
apnet2_ies_predicted = ap2.predict_qcel_mols(
    mols=df_plot['qcel_molecule'].tolist(),
)

# |%%--%%| <8jtfD3m3S4|nxmLPrfAfx>

# AP-Net2 IE
df_plot['APNet2'] = np.sum(apnet2_ies_predicted, axis=1)
df_plot['APNet2 error'] = (df_plot['APNet2'] - df_plot['reference']).astype(float)
print(df_plot[['APNet2', 'reference', 'PBE/6-31G*']])
print(apnet2_ies_predicted)
error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "B3LYP/6-31G*": "B3LYP/6-31G* error",
        "SAPT0/cc-pvdz": "SAPT0/cc-pvdz error",
        "APNet2": "APNet2 error",
    },
    output_filename="S22-IE-AP2.png",
    figure_size=(4, 4),
)

# |%%--%%| <nxmLPrfAfx|XtMJQEokjd>

# Training models on new QM data for Transfer Learning Task

from apnet_pt import pairwise_datasets

print(
    df_plot['qcel_molecule'].tolist(),
    df_plot['PBE/6-31G*'].tolist(),
)
print(len(df_plot['qcel_molecule'].tolist()))

ds2 = pairwise_datasets.apnet2_module_dataset(
    root="data_dir",
    spec_type=None,
    atom_model=apnet_pt.AtomModels.ap2_atom_model.AtomModel().set_pretrained_model(model_id=0),
    qcel_molecules=df_plot['qcel_molecule'].tolist(),
    energy_labels=[np.array([i]) for i in df_plot['reference'].tolist()],
    skip_compile=True,
    force_reprocess=True,
    atomic_batch_size=8,
    prebatched=False,
    in_memory=True,
    batch_size=2,
    # datapoint_storage_n_objects=4,
)
print("APNet2 dataset created")
print(ds2)
print(ds2.training_batch_size)
# print(ds2.data)

# |%%--%%| <XtMJQEokjd|xXslqNQSRI>

# Transfer Learning APNet2 model on computed QM data
print(ds2)
ap2.train(
    dataset=ds2,
    n_epochs=10,
    transfer_learning=True,
    skip_compile=True,
    model_path="apnet2_transfer_learning.pt",
    split_percent=0.8,
)
print(ap2)

# |%%--%%| <xXslqNQSRI|KPOeFMBhWm>

ap2_best_test = APNet2Model(
    atom_model=ap2.atom_model,
    pre_trained_model_path="apnet2_transfer_learning.pt",
)
print(ap2_best_test)

# AP-Net2 IE
apnet2_ies_predicted_final = ap2.predict_qcel_mols(
    mols=df_plot['qcel_molecule'].tolist(),
    batch_size=2,
)
apnet2_ies_predicted_best_test = ap2_best_test.predict_qcel_mols(
    mols=df_plot['qcel_molecule'].tolist(),
    batch_size=2,
)
print(apnet2_ies_predicted_final)
print(apnet2_ies_predicted_best_test)
df_plot['APNet2 final'] = np.sum(apnet2_ies_predicted_final, axis=1)
df_plot['APNet2 best test'] = np.sum(apnet2_ies_predicted_best_test, axis=1)
print(df_plot[['APNet2 final', 'APNet2 best test', 'reference', 'PBE/6-31G*']])
df_plot['APNet2 final error'] = (df_plot['APNet2 final'] - df_plot['reference']).astype(float)
df_plot['APNet2 best test error'] = (df_plot['APNet2 best test'] - df_plot['reference']).astype(float)

error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "B3LYP/6-31G*": "B3LYP/6-31G* error",
        "APNet2": "APNet2 error",
        "APNet2 final": "APNet2 final error",
        "APNet2 best test": "APNet2 best test error",
    },
    output_filename="S22-IE-AP2.png",
    figure_size=(4, 4),
)

# |%%--%%| <KPOeFMBhWm|OVVcYRXWUA>

# Be careful with this for it can corrupt running status...
!ps aux | grep qcfractal | awk '{ print $2 }' | xargs kill -9

# |%%--%%| <OVVcYRXWUA|ULuCE0r8ED>

# Load in a dataset from a recent Sherrill work (Levels of SAPT II)
df_LoS = pd.read_pickle("./combined_df_subset_358.pkl")
print(df_LoS[['Benchmark', 'SAPT2+3(CCD)DMP2 TOTAL ENERGY aqz', 'MP2 IE atz', 'SAPT0 TOTAL ENERGY adz' ]])

# Limit to 100 molecules with maximum of 16 atoms to keep computational cost down
df_LoS['size'] = df_LoS['atomic_numbers'].apply(lambda x: len(x))
df_LoS = df_LoS[df_LoS['size'] <= 16]
df_LoS = df_LoS.sample(100, random_state=42, axis=0).copy()
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

# |%%--%%| <ULuCE0r8ED|5HSrFMKRJh>


