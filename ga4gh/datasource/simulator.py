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
import ga4gh.protocol as protocol


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

    def run_get_bases(self, start, end):
        return self.bases[start:end]


class SimulatedVariantSet(registry.VariantSet):

    __tablename__ = 'SimulatedVariantSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('VariantSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity':'SimulatedVariantSet',
    }

    def __init__(
            self, name, random_seed, num_calls=1):
        super(SimulatedVariantSet, self).__init__(name)
        self.random_seed = random_seed
        self._create_metadata()
        for i in range(num_calls):
            name = "simCallSet_{}".format(i)
            call_set = registry.CallSet(name=name, sample_id=name)
            self.call_sets.append(call_set)

            # # build up infos of increasing size
            # for j in range(i):
            #     callSet._info["key_{}".format(j)] = "value_{}".format(j)

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

    def get_variant_id(self, variant):
        """
        Returns an ID suitable for the specifid variant.
        """
        # TODO this should be done generally using the CompoundId
        # infrastructure.
        return "{}:{}:{}".format(
            self.id, variant.reference_name, variant.start)

    def generate_variant(self, reference_name, call_sets, position):
        """
        Generate a random variant for the specified position.
        """
        rng = random.Random()
        rng.seed(self.random_seed + position)
        variant = protocol.Variant()
        variant.variant_set_id = str(self.id)
        variant.reference_name = reference_name
        variant.created = protocol.datetime_to_milliseconds(
            self.creation_timestamp)
        variant.updated = protocol.datetime_to_milliseconds(
            self.update_timestamp)
        variant.start = position
        variant.end = position + 1  # SNPs only for now
        bases = ["A", "C", "G", "T"]
        ref = rng.choice(bases)
        variant.reference_bases = ref
        alt = rng.choice([base for base in bases if base != ref])
        variant.alternate_bases.append(alt)
        for call_set in call_sets:
            call = variant.calls.add()
            call.call_set_id = str(call_set.id)
            call.call_set_name = call_set.name
            # for now, the genotype is either [0,1], [1,1] or [1,0] with equal
            # probability; probably will want to do something more
            # sophisticated later.
            random_choice = rng.choice([[0, 1], [1, 0], [1, 1]])
            call.genotype.extend(random_choice)
            # TODO What is a reasonable model for generating these likelihoods?
            # Are these log-scaled? Spec does not say.
            call.genotype_likelihood.extend([-100, -100, -100])
        variant.id = self.get_variant_id(variant)
        return variant

    def run_search(self, request, call_sets, response_builder):
        start = request.start
        if request.page_token:
            # TODO this should raise the correct error.
            start = int(request.page_token)
        i = start
        while i < request.end and not response_builder.isFull():
            response_builder.addValue(
                self.generate_variant(request.reference_name, call_sets, i))
            i += 1
        if i != request.end:
            response_builder.setNextPageToken(str(i))


#     def getVariant(self, compoundId):
#         randomNumberGenerator = random.Random()
#         start = int(compoundId.start)
#         randomNumberGenerator.seed(self._randomSeed + start)
#         variant = self.generateVariant(
#             compoundId.reference_name, start, randomNumberGenerator)
#         return variant

#     def getVariants(self, referenceName, startPosition, endPosition,
#                     callSetIds=None):
#         randomNumberGenerator = random.Random()
#         randomNumberGenerator.seed(self._randomSeed)
#         i = startPosition
#         while i < endPosition:
#             if randomNumberGenerator.random() < self._variantDensity:
#                 randomNumberGenerator.seed(self._randomSeed + i)
#                 yield self.generateVariant(
#                     referenceName, i, randomNumberGenerator)
#             i += 1




