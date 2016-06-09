"""
Module defining concrete instances of the top level containers that generate
simulated data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import hashlib

import sqlalchemy
import sqlalchemy.orm as orm

import ga4gh.registry as registry


class SimulatedReferenceSet(registry.ReferenceSet):
    __tablename__ = 'SimulatedReferenceSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReferenceSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity':'SimulatedReferenceSet',
    }

    def __init__(
            self, name, random_seed, num_references=10, reference_length=200):
        super(SimulatedReferenceSet, self).__init__(name)
        self.random_seed = random_seed
        rng = random.Random()
        rng.seed(random_seed)

        self.description = "Simulated reference set"
        self.assemblyId = str(random.randint(0, 2**32))
        self.isDerived = bool(random.randint(0, 1))
        self.ncbiTaxonId = random.randint(0, 2**16)
        # self._sourceAccessions = []
        # for i in range(random.randint(1, 3)):
        #         self._sourceAccessions.append("sim_accession_{}".format(
        #             random.randint(1, 2**32)))
        self._sourceUri = "http://example.com/reference.fa"
        for i in range(num_references):
            reference_seed = rng.getrandbits(32)
            reference_name = "srs{}".format(i)
            reference = SimulatedReference(
                reference_name, reference_seed, reference_length)
            self.references.append(reference)
        self._compute_md5checksum()

class SimulatedReference(registry.Reference):
    __tablename__ = 'SimulatedReference'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('Reference.id'),
        primary_key=True)
    bases = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity':'SimulatedReference',
    }

    def __init__(self, name, random_seed, length):
        super(SimulatedReference, self).__init__(name)
        rng = random.Random()
        rng.seed(random_seed)
        self.length = length
        self.bases = "".join([rng.choice('ACGT') for _ in range(length)])
        self.md5checksum = hashlib.md5(self.bases).hexdigest()
        self.is_derived = bool(rng.randint(0, 1))
        self.source_divergence = 0
        if self.is_derived:
            self.source_divergence = rng.uniform(0, 0.1)
        # self._sourceAccessions = []
        # for i in range(random.randint(1, 3)):
        #         self._sourceAccessions.append("sim_accession_{}".format(
        #             random.randint(1, 2**32)))
        self._sourceUri = "http://example.com/reference.fa"

    def getBases(self, start, end):
        self.checkQueryRange(start, end)
        return self._bases[start:end]


class SimulatedVariantSet(registry.VariantSet):

    __tablename__ = 'SimulatedVariantSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('VariantSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    variant_density = sqlalchemy.Column(sqlalchemy.Float, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity':'SimulatedVariantSet',
    }

    def __init__(
            self, name, random_seed, num_calls=1, variant_density=1):
        super(SimulatedVariantSet, self).__init__(name)
        self.random_seed = random_seed
        self.variant_density = variant_density

        self._create_metadata()
        print("New simulated variant set")


    def _create_metadata(self):
        metadata = registry.VariantSetMetadata(
            key="version", value="VCFv4.1", type="String", number="1")
        self.variant_set_metadata.append(metadata)
        metadata = registry.VariantSetMetadata(
            key="INFO.FIELD1", value="", type="String",
            description="FIELD1 description")
        self.variant_set_metadata.append(metadata)
        metadata = registry.VariantSetMetadata(
            key="INFO.FIELD2", type="Integer", number="1",
            description="FIELD2 description")
        self.variant_set_metadata.append(metadata)


    def run_search(self, request, response_builder):
        print("running search on ", request, response_builder)
        for j in range(10):
            variant = protocol.Variant()
            variant.variant_set_id = str(self.id)
            variant.reference_name = request.reference_name
            variant.start = j
            variant.end = j + 1
            response_builder.addValue(variant)




