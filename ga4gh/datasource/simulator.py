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

import ga4gh.registry as registry
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


#####################################################################
#
# References
#
#####################################################################


class SimulatedReferenceSet(registry.ReferenceSet):
    """
    A simulated reference set. Generates values for class attributes randomly
    and creates a given number of SimulatedReference instances.
    """
    __tablename__ = 'SimulatedReferenceSet'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReferenceSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'SimulatedReferenceSet',
    }

    def __init__(
            self, name, random_seed, num_references=10, reference_length=200):
        super(SimulatedReferenceSet, self).__init__(name)
        self.random_seed = random_seed
        rng = random.Random()
        rng.seed(random_seed)

        self.description = "Simulated reference set"
        self.assembly_id = str(random.randint(0, 2**32))
        self.is_derived = bool(random.randint(0, 1))
        self.ncbi_taxon_id = random.randint(0, 2**16)
        for i in range(random.randint(1, 3)):
            accession = registry.Accession(
                "sim_accession_{}".format(random.randint(1, 2**32)))
            self.source_accessions.append(accession)
        self._sourceUri = "http://example.com/reference.fa"
        for i in range(num_references):
            reference_seed = rng.getrandbits(32)
            reference_name = "srs{}".format(i)
            reference = SimulatedReference(
                reference_name, reference_seed, reference_length)
            self.references.append(reference)
        self._compute_md5checksum()


class SimulatedReference(registry.Reference):
    """
    A simulated reference. Stores a random sequence of a given length, and
    generates remaining attributes randomly.
    """
    __tablename__ = 'SimulatedReference'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('Reference.id'),
        primary_key=True)
    bases = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'SimulatedReference',
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
        for i in range(random.randint(1, 3)):
            accession = registry.Accession(
                "sim_accession_{}".format(random.randint(1, 2**32)))
            self.source_accessions.append(accession)
        self._sourceUri = "http://example.com/reference.fa"

    def run_get_bases(self, start, end):
        return self.bases[start:end]


#####################################################################
#
# Variants
#
#####################################################################


class SimulatedVariantSet(registry.VariantSet):

    __tablename__ = 'SimulatedVariantSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('VariantSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'SimulatedVariantSet',
    }

    def __init__(
            self, name, random_seed, bio_samples):
        super(SimulatedVariantSet, self).__init__(name)
        self.random_seed = random_seed
        self._create_metadata()
        for i, bio_sample in enumerate(bio_samples):
            name = "simCallSet_{}".format(i)
            call_set = registry.CallSet(name=name)
            call_set.bio_sample = bio_sample
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
            try:
                start = int(request.page_token)
            except ValueError:
                raise exceptions.BadPageTokenException(request.page_token)
        i = start
        while i < request.end and not response_builder.isFull():
            response_builder.addValue(
                self.generate_variant(request.reference_name, call_sets, i))
            i += 1
        if i != request.end:
            response_builder.setNextPageToken(str(i))

    def get_variants(self, reference_name, start=0, end=None, call_sets=[]):
        """
        Returns an iterator over the variants in the specified region.
        """
        the_end = end
        if end is None:
            the_end = 2**32
        for i in range(start, the_end):
            yield self.generate_variant(reference_name, call_sets, i)


#####################################################################
#
# Reads
#
#####################################################################


class SimulatedReadGroupSet(registry.ReadGroupSet):

    __tablename__ = 'SimulatedReadGroupSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReadGroupSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'SimulatedReadGroupSet',
    }

    def __init__(self, name, bio_samples, random_seed=1, num_read_groups=2):
        super(SimulatedReadGroupSet, self).__init__(name)
        self.random_seed = random_seed
        rng = random.Random()
        rng.seed(random_seed)
        for j in range(num_read_groups):
            read_group = SimulatedReadGroup(
                name="sim_rgs{}".format(j),
                random_seed=rng.randint(0, 2**31))
            read_group.bio_sample = rng.choice(bio_samples)
            read_group.sample_name = read_group.bio_sample.individual.name
            self.read_groups.append(read_group)

    def generate_alignment(self, reference, read_groups, start):
        # TODO fill out a bit more
        rng = random.Random(self.random_seed + start)
        alignment = protocol.ReadAlignment()
        alignment.read_group_id = str(rng.choice(read_groups).id)
        alignment.fragment_length = rng.randint(10, 100)
        alignment.aligned_sequence = ""
        for i in range(alignment.fragment_length):
            # TODO: are these reasonable quality values?
            alignment.aligned_quality.append(rng.randint(1, 20))
            alignment.aligned_sequence += rng.choice("ACGT")

        alignment.alignment.position.position = 0
        alignment.alignment.position.reference_name = "NotImplemented"
        alignment.alignment.position.strand = protocol.POS_STRAND
        alignment.duplicate_fragment = False
        alignment.failed_vendor_quality_checks = False

        alignment.fragment_name = "{}$simulated{}".format(self.name, i)
        alignment.number_reads = 0
        alignment.improper_placement = False
        alignment.read_number = -1
        alignment.secondary_alignment = False
        alignment.supplementary_alignment = False
        # alignment.id = self._parentContainer.getReadAlignmentId(alignment)
        return alignment

    def run_search(self, request, reference, read_groups, response_builder):
        # TODO abstract this so that we can share this paging code with
        # the variant set.
        start = request.start
        if request.page_token:
            try:
                start = int(request.page_token)
            except ValueError:
                raise exceptions.BadPageTokenException(request.page_token)
        i = start
        while i < request.end and not response_builder.isFull():
            response_builder.addValue(
                self.generate_alignment(reference, read_groups, i))
            i += 1
        if i != request.end:
            response_builder.setNextPageToken(str(i))


class SimulatedReadGroup(registry.ReadGroup):

    __tablename__ = 'SimulatedReadGroup'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReadGroup.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'SimulatedReadGroup',
    }

    def __init__(self, name, random_seed=1):
        super(SimulatedReadGroup, self).__init__(name=name)
        self.random_seed = random_seed
        rng = random.Random()
        rng.seed(random_seed)
        self.predicted_insert_size = random.randint(10, 100)
        self.read_stats = registry.ReadStats(
            aligned_read_count=random.randint(0, 100),
            unaligned_read_count=random.randint(0, 100),
            base_count=random.randint(0, 100))

#####################################################################
#
# Features.
#
#####################################################################


class SimulatedFeatureSet(registry.FeatureSet):
    """
    Simulated data backend for FeatureSet, used for internal testing.
    """

    __tablename__ = 'SimulatedFeatureSet'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('FeatureSet.id'),
        primary_key=True)
    random_seed = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity': 'SimulatedFeatureSet',
    }

    def __init__(self, name, random_seed=1):
        super(SimulatedFeatureSet, self).__init__(name)
        self.random_seed = random_seed

    def _getRandomfeatureType(self, randomNumberGenerator):
        ontologyTuples = [
            ("gene", "SO:0000704"),
            ("exon", "SO:0000147")]
        term = protocol.OntologyTerm()
        ontologyTuple = randomNumberGenerator.choice(ontologyTuples)
        term.term, term.id = ontologyTuple[0], ontologyTuple[1]
        term.source_name = "sequenceOntology"
        term.source_version = "0"
        return term

    def _generateSimulatedFeature(self, randomNumberGenerator):
        feature = protocol.Feature()
        feature.feature_set_id = self.getId()
        feature.start = randomNumberGenerator.randint(1000, 2000)
        feature.end = feature.start + randomNumberGenerator.randint(1, 100)
        feature.feature_type.CopyFrom(self._getRandomfeatureType(
            randomNumberGenerator))
        strands = [protocol.POS_STRAND, protocol.NEG_STRAND]
        feature.strand = randomNumberGenerator.choice(strands)
        attributes = {
            "gene_name": "Frances",
            "gene_type": "mRNA",
            "gene_status": "UNKNOWN"}
        for key, value in attributes.items():
            feature.attributes.vals[key].values.add().string_value = value
        return feature

    def getFeature(self, compoundId):
        """
        Fetches a simulated feature by ID.

        :param compoundId: any non-null string
        :return: A simulated feature with id set to the same value as the
            passed-in compoundId.
        ":raises: exceptions.ObjectWithIdNotFoundException if None is passed
            in for the compoundId.
        """
        if compoundId is None:
            raise exceptions.ObjectWithIdNotFoundException(compoundId)
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        feature = self._generateSimulatedFeature(randomNumberGenerator)
        feature.id = str(compoundId)
        feature.parent_id = ""  # TODO: Test with nonempty parentIDs?
        return feature

    def getFeatures(
            self, referenceName, start, end,
            pageToken, pageSize,
            featureTypes=[], parentId=None, numFeatures=10):
        """
        Returns a set number of simulated features.

        :param referenceName: name of reference to "search" on
        :param start: start coordinate of query
        :param end: end coordinate of query
        :param pageToken: None or int
        :param pageSize: None or int
        :param featureTypes: optional list of ontology terms to limit query
        :param parentId: optional parentId to limit query.
        :param numFeatures: number of features to generate in the return.
            10 is a reasonable (if arbitrary) default.
        :return: Yields feature, nextPageToken pairs.
            nextPageToken is None if last feature was yielded.
        """
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        if pageToken:
            nextPageToken = int(pageToken)
        else:
            nextPageToken = 0
        for featureId in range(numFeatures):
            gaFeature = self._generateSimulatedFeature(randomNumberGenerator)
            gaFeature.id = self.getCompoundIdForFeatureId(featureId)
            match = (
                gaFeature.start < end and
                gaFeature.end > start and
                gaFeature.reference_name == referenceName and (
                    featureTypes is None or len(featureTypes) == 0 or
                    gaFeature.feature_type in featureTypes))
            if match:
                gaFeature.parent_id = ""  # TODO: Test nonempty parentIDs?
                if nextPageToken < numFeatures - 1:
                    nextPageToken += 1
                else:
                    nextPageToken = None
                yield gaFeature, (
                    str(nextPageToken)
                    if nextPageToken is not None else None)

    def generate_feature(self, reference_name, position):
        rng = random.Random()
        rng.seed(self.random_seed + position)
        feature = protocol.Feature()
        feature.feature_set_id = str(self.id)
        feature.start = position
        feature.end = feature.start + rng.randint(1, 100)
        feature.feature_type.CopyFrom(self._getRandomfeatureType(rng))
        feature.reference_name = reference_name
        strands = [protocol.POS_STRAND, protocol.NEG_STRAND]
        feature.strand = rng.choice(strands)
        attributes = {
            "gene_name": "Frances",
            "gene_type": "mRNA",
            "gene_status": "UNKNOWN"}
        for key, value in attributes.items():
            feature.attributes.vals[key].values.add().string_value = value
        return feature

    def run_search(self, request, response_builder):
        # TODO abstract this code so it can be shared with other simulators.
        start = request.start
        if request.page_token:
            try:
                start = int(request.page_token)
            except ValueError:
                raise exceptions.BadPageTokenException(request.page_token)
        i = start
        while i < request.end and not response_builder.isFull():
            response_builder.addValue(
                self.generate_feature(request.reference_name, i))
            i += 1
        if i != request.end:
            response_builder.setNextPageToken(str(i))
