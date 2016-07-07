from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob

import ga4gh.registry as registry
import ga4gh.datasource.simulator as simulator
import ga4gh.datasource.htslib as htslib

db = registry.RegistryDb("postgres://ga4gh-dev:password@localhost/ga4gh-registry")
# db = registry.RegistryDb("sqlite:///registry.db")
db.open()
db.initialise()

rs1 = htslib.HtslibReferenceSet(
    "NCBI37", "tests/data/referenceSets/NCBI37.fa.gz")
rs1.ncbi_taxon_id = 9606
rs1.is_derived = False
db.add_reference_set(rs1)

for j in range(1):
    rs2 = simulator.SimulatedReferenceSet("simrs_{}".format(j), j)
    db.add_reference_set(rs2)

# db.commit()

dataset = registry.Dataset("ds1")
db.add_dataset(dataset)

vcfs = glob.glob(
    "tests/data/datasets/dataset1/variants/1kgPhase1/chr*.vcf.gz")
indexes = [vcf + ".tbi" for vcf in vcfs]

variant_set = htslib.HtslibVariantSet("vs1", vcfs, indexes)
variant_set.dataset = dataset
variant_set.reference_set = rs1
db.add_variant_set(variant_set)

variant_set = simulator.SimulatedVariantSet("sim_vs", 234, num_calls=10)
variant_set.dataset = dataset
variant_set.reference_set = rs1
db.add_variant_set(variant_set)

print(variant_set.get_protobuf())

for metadata in variant_set.variant_set_metadata:
    print(metadata.id, metadata.key)

for call_set in variant_set.call_sets:
    print(call_set.name, call_set.variant_sets)
    print(call_set.get_protobuf())

for rs in db.get_reference_sets():
    print(rs.id, rs.name)

# rgs = registry.ReadGroupSet("name")
# rgs.dataset = dataset
# rgs.reference_set = rs1
# rgs.read_stats = registry.ReadStats(1, 13, 13)
# rg = registry.ReadGroup("rg1", "sid")
# rgs.read_groups.append(rg)
# rg.read_stats = registry.ReadStats(1, 1, 1)
# program = registry.Program(name="x", version="y")
# rgs.programs.append(program)
# rgs.experiment = registry.Experiment(description="fake")
# rgs = simulator.SimulatedReadGroupSet(
#     "sim_rgs", 100, num_read_groups=5)
# rgs.dataset = dataset
# rgs.reference_set = rs1

# db.add_read_group_set(rgs)

# bam = "tests/data/datasets/dataset1/reads/wgEncodeUwRepliSeqBg02esG1bAlnRep1_sample.bam"
# rgs = htslib.HtslibReadGroupSet("wg", bam, bam + ".bai")
# rgs.dataset = dataset
# rgs.reference_set = rs1
# db.add_read_group_set(rgs)

bam = "tests/data/datasets/dataset1/reads/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
rgs = htslib.HtslibReadGroupSet("HG0096", bam, bam + ".bai")
rgs.dataset = dataset
rgs.reference_set = rs1
db.add_read_group_set(rgs)

# bam = "1kg-data/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
# rgs = htslib.HtslibReadGroupSet("NA12878", bam, bam + ".bai")
# rgs.dataset = dataset
# rgs.reference_set = rs1
# db.add_read_group_set(rgs)

# bam = "1kg-data/NA21144.mapped.ILLUMINA.bwa.GIH.low_coverage.20130415.bam"
# rgs = htslib.HtslibReadGroupSet("NA21144", bam, bam + ".bai")
# rgs.dataset = dataset
# rgs.reference_set = rs1
# db.add_read_group_set(rgs)



# print(rgs)
# # print(rgs.dataset)
# print(rgs.get_protobuf())

db.close()
