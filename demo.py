
import psi4
from qcportal import PortalClient
from pprint import pprint as pp
from qcelemental.models import Molecule
import qcelemental as qcel
from qcportal.singlepoint import SinglepointDataset, SinglepointDatasetEntry, QCSpecification
import pandas as pd
import numpy as np

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
        "service_frequency": 10,
        "max_active_services": None,
        "heartbeat_frequency": None,
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
for rec in client.query_records():
    pp(rec)

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

for rec in client.query_records():
    pp(rec.dict)
    pp(rec.error)

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

# |%%--%%| <dSW1A9HxYB|dwYb9dbQNI>

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
    charges = qcvars["mbis charges"]
    dipoles = qcvars["mbis dipoles"]
    quadrupoles = qcvars["mbis quadrupoles"]
    level_of_theory = f"{record_dict['specification']['method']}/{record_dict['specification']['basis']}"

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
        "molecular_multiplicity"
    ]

df = ds.compile_values(
    value_call=assemble_multipole_data,
    value_names=assemble_multipole_data_value_names(),
    unpack=True,
)
print(df)
pp(df.columns.tolist())
print(df['psi4/hf/sto-3g'])
pp(df['psi4/hf/sto-3g'].columns.tolist())
df_hf_sto3g = df['psi4/hf/sto-3g']
print(df_hf_sto3g.columns.tolist())
print(df_hf_sto3g)

#|%%--%%| <7RHL31QOoC|CDD1QjxHpc>

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

#|%%--%%| <CDD1QjxHpc|XT87RegBfm>

from cdsg_plot import error_statistics

pp(df_mb.columns.tolist())

df_plot = pd.DataFrame(
    {
        "HF/6-31G*": df_mb["psi4/hf/6-31g*"]["CP_IE"],
        "PBE/6-31G*": df_mb["psi4/pbe/6-31g*"]["CP_IE"],
        "B3LYP/6-31G*": df_mb["psi4/b3lyp/6-31g*"]["CP_IE"],
        'reference': ref_IEs,
    }
)
df_plot['HF/6-31G* error'] = df_plot['HF/6-31G*'] - df_plot['reference']
df_plot['PBE/6-31G* error'] = df_plot['PBE/6-31G*'] - df_plot['reference']
df_plot['B3LYP/6-31G* error'] = df_plot['B3LYP/6-31G*'] - df_plot['reference']
pd.set_option('display.max_rows', None)
print(df_plot)
print(df_plot.describe())


#|%%--%%| <XT87RegBfm|dRiuyCtOh1>


error_statistics.violin_plot(
    df_plot,
    df_labels_and_columns={
        "HF/6-31G*": "HF/6-31G* error",
        "PBE/6-31G*": "PBE/6-31G* error",
        "B3LYP/6-31G*": "B3LYP/6-31G* error",
    },
    output_filename="S22-IE.png",
)

#|%%--%%| <dRiuyCtOh1|fLqWKAJRoW>
r"""°°°
![S22-IE_violin.png](./S22-IE_violin.png)
°°°"""
# |%%--%%| <fLqWKAJRoW|OVVcYRXWUA>

# Be careful with this for it can corrupt running status...
# !ps aux | grep qcfractal | awk '{ print $2 }' | xargs kill -9

# |%%--%%| <OVVcYRXWUA|5HSrFMKRJh>


