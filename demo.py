
import psi4
from qcportal import PortalClient
from pprint import pprint as pp
from qcelemental.models import Molecule
from qcportal.singlepoint import SinglepointDataset, SinglepointDatasetEntry, QCSpecification

# need manybodydataset
from qcportal.manybody import ManybodyDataset, ManybodyDatasetEntry, ManybodyDatasetSpecification

#|%%--%%| <virdG1SNXY|BVc6W6uOta>

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
)


#|%%--%%| <BVc6W6uOta|CwEhpqwLXX>

!qcfractal-server --config=`pwd`/qcfractal/qcfractal_config.yaml start > qcfractal_server.log & disown

# NOTE kill server when finished by running:
#     ps aux | grep qcfractal-server | awk '{ print $2 }'
#     kill -9 <PID>

#|%%--%%| <CwEhpqwLXX|3HjtiyIuFg>

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
#|%%--%%| <hMCmRgdGJ4|Z0wXrcgRq8>

# Now create S22 Interaction Energy Dataset
from s22 import geoms

# geoms is a list of qcelemental Molecule objects that can be used to create a
# QCArchive dataset
print(geoms)

#|%%--%%| <Z0wXrcgRq8|i8ICwzPWaD>

# Create client dataset

ds_name = 'S22-multipoles'

try:
    ds = client.add_dataset("singlepoint", ds_name,
                            f"Dataset to contain {ds_name}")
    print(f"Added {ds_name} as dataset")
except Exception:
    ds = client.get_dataset("singlepoint", ds_name)
    print(f"Found {ds_name} dataset, using this instead")
    print(ds)

#|%%--%%| <i8ICwzPWaD|1mV2xNrvP7>

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

#|%%--%%| <1mV2xNrvP7|TdmQndJlQR>

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

#|%%--%%| <TdmQndJlQR|eifKS4ffOR>

# Run the computations
ds.submit()
print(f"Submitted {ds_name} dataset")

#|%%--%%| <eifKS4ffOR|vsEEQZOLgv>

# Check the status of the dataset - can repeatedly run this to see the progress
ds.status()

#|%%--%%| <vsEEQZOLgv|g31JlHrgso>



#|%%--%%| <g31JlHrgso|OVVcYRXWUA>

!ps aux | grep qcfractal-server | awk '{ print $2 }'
#|%%--%%| <OVVcYRXWUA|5HSrFMKRJh>



