"""
Module defining concrete instances of the top level containers using
standard data files interpreted using pysam/htslib.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import hashlib
import os

import pysam
import sqlalchemy
import sqlalchemy.orm as orm
import google.protobuf.struct_pb2 as struct_pb2

import ga4gh.registry as registry
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


DEFAULT_READGROUP_NAME = "Default_RG"
"""
This is the name used for the default read group created when there
are no RG headers in a BAM file.
"""


_nothing = object()


def is_empty_iter(it):
    """
    Return True iff the iterator is empty or exhausted
    """
    return next(it, _nothing) is _nothing


def _encodeValue(value):
    # TODO give this a better name.
    if isinstance(value, (list, tuple)):
        return [struct_pb2.Value(string_value=str(v)) for v in value]
    else:
        return [struct_pb2.Value(string_value=str(value))]


def parseMalformedBamHeader(headerDict):
    """
    Parses the (probably) intended values out of the specified
    BAM header dictionary, which is incompletely parsed by pysam.
    This is caused by some tools incorrectly using spaces instead
    of tabs as a seperator.
    """
    headerString = " ".join(
        "{}:{}".format(k, v) for k, v in headerDict.items() if k != 'CL')
    ret = {}
    for item in headerString.split():
        key, value = item.split(":", 1)
        # build up dict, casting everything back to original type
        ret[key] = type(headerDict.get(key, ""))(value)
    if 'CL' in headerDict:
        ret['CL'] = headerDict['CL']
    return ret


class SamCigar(object):
    """
    Utility class for working with SAM CIGAR strings
    """
    # see http://pysam.readthedocs.org/en/latest/api.html
    # #pysam.AlignedSegment.cigartuples
    cigarStrings = [
        protocol.CigarUnit.ALIGNMENT_MATCH,
        protocol.CigarUnit.INSERT,
        protocol.CigarUnit.DELETE,
        protocol.CigarUnit.SKIP,
        protocol.CigarUnit.CLIP_SOFT,
        protocol.CigarUnit.CLIP_HARD,
        protocol.CigarUnit.PAD,
        protocol.CigarUnit.SEQUENCE_MATCH,
        protocol.CigarUnit.SEQUENCE_MISMATCH,
    ]

    @classmethod
    def ga2int(cls, value):
        for i, cigarString in enumerate(cls.cigarStrings):
            if value == cigarString:
                return i

    @classmethod
    def int2ga(cls, value):
        return cls.cigarStrings[value]


class SamFlags(object):
    """
    Utility class for working with SAM flags
    """
    READ_PAIRED = 0x1
    READ_PROPER_PAIR = 0x2
    READ_UNMAPPED = 0x4
    MATE_UNMAPPED = 0x8
    READ_REVERSE_STRAND = 0x10
    MATE_REVERSE_STRAND = 0x20
    FIRST_IN_PAIR = 0x40
    SECOND_IN_PAIR = 0x80
    SECONDARY_ALIGNMENT = 0x100
    FAILED_QUALITY_CHECK = 0x200
    DUPLICATE_READ = 0x400
    SUPPLEMENTARY_ALIGNMENT = 0x800

    @staticmethod
    def isFlagSet(flagAttr, flag):
        return flagAttr & flag == flag

    @staticmethod
    def setFlag(flagAttr, flag):
        return flagAttr | flag


class PysamMixin(object):
    """
    A mixin class to simplify working with DatamodelObjects based on
    directories of files interpreted using pysam. This mixin is designed
    to work within the DatamodelObject hierarchy.
    """
    samMin = 0
    samMaxStart = 2**30 - 1
    samMaxEnd = 2**30

    vcfMin = -2**31
    vcfMax = 2**31 - 1

    fastaMin = 0
    fastaMax = 2**30 - 1

    rNameMin = 0
    rNameMax = 85

    maxStringLength = 2**10  # arbitrary

    @classmethod
    def sanitizeVariantFileFetch(cls, contig=None, start=None, stop=None):
        if contig is not None:
            contig = cls.sanitizeString(contig, 'contig')
        if start is not None:
            start = cls.sanitizeInt(start, cls.vcfMin, cls.vcfMax, 'start')
        if stop is not None:
            stop = cls.sanitizeInt(stop, cls.vcfMin, cls.vcfMax, 'stop')
        if start is not None and stop is not None:
            cls.assertValidRange(start, stop, 'start', 'stop')
        return contig, start, stop

    @classmethod
    def sanitizeAlignmentFileFetch(cls, start=None, end=None):
        if start is not None:
            start = cls.sanitizeInt(
                start, cls.samMin, cls.samMaxStart, 'start')
        if end is not None:
            end = cls.sanitizeInt(end, cls.samMin, cls.samMaxEnd, 'end')
        if start is not None and end is not None:
            cls.assertValidRange(start, end, 'start', 'end')
        return start, end

    @classmethod
    def assertValidRange(cls, start, end, startName, endName):
        if start > end:
            message = "invalid coordinates: {} ({}) " \
                "greater than {} ({})".format(startName, start, endName, end)
            raise exceptions.DatamodelValidationException(message)

    @classmethod
    def assertInRange(cls, attr, minVal, maxVal, attrName):
        message = "invalid {} '{}' outside of range [{}, {}]"
        if attr < minVal:
            raise exceptions.DatamodelValidationException(message.format(
                attrName, attr, minVal, maxVal))
        if attr > maxVal:
            raise exceptions.DatamodelValidationException(message.format(
                attrName, attr, minVal, maxVal))

    @classmethod
    def assertInt(cls, attr, attrName):
        if not isinstance(attr, (int, long)):
            message = "invalid {} '{}' not an int".format(attrName, attr)
            raise exceptions.DatamodelValidationException(message)

    @classmethod
    def sanitizeInt(cls, attr, minVal, maxVal, attrName):
        cls.assertInt(attr, attrName)
        if attr < minVal:
            attr = minVal
        if attr > maxVal:
            attr = maxVal
        return attr

    @classmethod
    def sanitizeString(cls, attr, attrName):
        if not isinstance(attr, basestring):
            message = "invalid {} '{}' not a string".format(
                attrName, attr)
            raise exceptions.DatamodelValidationException(message)
        if isinstance(attr, unicode):
            attr = attr.encode('utf8')
        if len(attr) > cls.maxStringLength:
            attr = attr[:cls.maxStringLength]
        return attr

    def get_file_handle(self, open_func, url, *args, **kwargs):
        return file_handle_cache.get_file_handle(
            url, open_func, url, *args, **kwargs)


class HtslibReference(registry.Reference, PysamMixin):
    __tablename__ = 'HtslibReference'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('Reference.id'),
        primary_key=True)
    fasta_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    __mapper_args__ = {
        'polymorphic_identity': 'HtslibReference',
    }

    def run_get_bases(self, start, end):
        self.check_query_range(start, end)
        fasta_file = self.get_file_handle(pysam.FastaFile, self.fasta_url)
        bases = fasta_file.fetch(self.name.encode(), start, end)
        return bases


class HtslibReferenceSet(registry.ReferenceSet):
    __tablename__ = 'HtslibReferenceSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReferenceSet.id'),
        primary_key=True)

    __mapper_args__ = {
        'polymorphic_identity': 'HtslibReferenceSet',
    }

    def __init__(self, name, fasta_url):
        super(HtslibReferenceSet, self).__init__(name)
        with pysam.FastaFile(fasta_url) as fasta_file:
            for reference_name in fasta_file.references:
                reference = HtslibReference(reference_name)
                bases = fasta_file.fetch(reference_name)
                reference.md5checksum = hashlib.md5(bases).hexdigest()
                reference.length = len(bases)
                reference.fasta_url = fasta_url
                self.references.append(reference)
        # Now that we have all the references, set the md5checksum
        self._compute_md5checksum()

##################################################
#
# Variants
#
##################################################


class HtslibVariantCompoundId(registry.VariantCompoundId):
    # TODO Implement this and document when we're fixing the CompoundId
    # system.
    pass


class VcfFile(registry.SqlAlchemyBase):
    __tablename__ = 'VcfFile'

    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    data_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    index_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    reference_name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    variant_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('HtslibVariantSet.id'))

    def __init__(self, data_url, index_url, reference_name):
        self.data_url = data_url
        self.index_url = index_url
        self.reference_name = reference_name


class HtslibVariantSet(registry.VariantSet, PysamMixin):
    __tablename__ = 'HtslibVariantSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('VariantSet.id'),
        primary_key=True)

    # Each VCF file row contains the reference_name to allow us to map
    # back to the file that contains a particular reference at query time.
    # Some VCF files contain large numbers of references, so we don't want
    # to load up this entire least on each query. Dynamic collections
    # allow us to explicitly query for the file we want.
    vcf_files = orm.relationship(
        "VcfFile", cascade="all, delete-orphan", lazy="dynamic")

    __mapper_args__ = {
        'polymorphic_identity': 'HtslibVariantSet',
    }

    def __init__(self, name, data_urls, index_urls):
        super(HtslibVariantSet, self).__init__(name)
        self.vcf_files = []
        for data_url, index_file in zip(data_urls, index_urls):
            var_file = pysam.VariantFile(data_url, index_filename=index_file)
            try:
                self._add_vcf_file(var_file, data_url, index_file)
            finally:
                var_file.close()
        self._check_state()

    def _check_state(self):
        """
        Checks the state of this variant set to ensure it's consistent before
        we commit it to the DB.
        """
        reference_names = set()
        for vcf_file in self.vcf_files:
            reference_name = vcf_file.reference_name
            if reference_name in reference_names:
                raise exceptions.OverlappingVcfException(
                    vcf_file.data_url, reference_name)
            reference_names.add(reference_name)

    def _add_vcf_file(self, var_file, data_url, index_file):
        """
        Updates the state of this variant set based on the information in the
        specified pysam VariantFile object.
        """
        if var_file.index is None:
            raise exceptions.NotIndexedException(data_url)
        for reference_name in var_file.index:
            # Unlike Tabix indices, CSI indices include all contigs defined
            # in the BCF header.  Thus we must test each one to see if
            # records exist or else they are likely to trigger spurious
            # overlapping errors.
            if not is_empty_iter(var_file.fetch(reference_name)):
                # We add an instance of VcfFile for each reference within the
                # VCF because we need to find the file that maps to a
                # particular reference.
                self.vcf_files.append(
                    VcfFile(data_url, index_file, reference_name))

        self._update_metadata(var_file)
        self._updateCallSetIds(var_file)
        # self._updateVariantAnnotationSets(varFile, dataUrl)

    def _update_metadata(self, var_file):
        if len(self.variant_set_metadata) == 0:
            self.variant_set_metadata = self._getMetadataFromVcf(var_file)

    def _updateCallSetIds(self, variantFile):
        """
        Updates the call set IDs based on the specified variant file.
        """
        if self.call_sets.count() == 0:
            for sample in variantFile.header.samples:
                call_set = registry.CallSet(name=sample, sample_id=sample)
                self.call_sets.append(call_set)
        # TODO check the callsets are already in here as required, or throw
        # an error.

    def _getMetadataFromVcf(self, varFile):
        # TODO document this function and refactor it to use snake_case and
        # get rid of the redundant buildMetadata function.

        # TODO replace buildMetadata with constructor for VariantSetMetadata
        def buildMetadata(
                key, type_="String", number="1", value="", description=""):
            metadata = registry.VariantSetMetadata()
            metadata.key = key
            metadata.value = value
            metadata.type = type_
            metadata.number = number
            metadata.description = description
            return metadata

        header = varFile.header
        ret = []
        ret.append(buildMetadata(key="version", value=header.version))
        formats = header.formats.items()
        infos = header.info.items()
        # TODO: currently ALT field is not implemented through pysam
        # NOTE: contigs field is different between vcf files,
        # so it's not included in metadata
        # NOTE: filters in not included in metadata unless needed
        for prefix, content in [("FORMAT", formats), ("INFO", infos)]:
            for contentKey, value in content:
                description = value.description.strip('"')
                key = "{0}.{1}".format(prefix, value.name)
                if key != "FORMAT.GT":
                    ret.append(buildMetadata(
                        key=key, type_=value.type,
                        number="{}".format(value.number),
                        description=description))
        return ret

    def _update_ga_call(self, ga_call, pysam_call):
        phaseset = ""
        if pysam_call.phased:
            # TODO this gives phaseset = "True"; is this correct???
            phaseset = str(pysam_call.phased)
        genotypeLikelihood = []
        info = {}
        for key, value in pysam_call.iteritems():
            if key == 'GL' and value is not None:
                genotypeLikelihood = list(value)
            elif key != 'GT':
                info[key] = _encodeValue(value)
        ga_call.genotype.extend(list(pysam_call.allele_indices))
        ga_call.phaseset = phaseset
        ga_call.genotype_likelihood.extend(genotypeLikelihood)
        for key in info:
            ga_call.info[key].values.extend(info[key])

    def get_variant_id(self, variant):
        bases = variant.reference_bases + str(tuple(variant.alternate_bases))
        bases_hash = hashlib.md5(bases).hexdigest()
        # TODO get rid of redundant None here
        compound_id = registry.VariantCompoundId(
            None, str(self.id), variant.reference_name, str(variant.start),
            bases_hash)
        return str(compound_id)

    def convert_variant(self, record, call_sets):
        """
        Converts the specified pysam variant record into a GA4GH Variant
        object. Only calls for the specified list of call_sets will
        be included.
        """
        variant = protocol.Variant()
        variant.variant_set_id = str(self.id)
        variant.reference_name = record.contig
        variant.created = protocol.datetime_to_milliseconds(
            self.creation_timestamp)
        variant.updated = protocol.datetime_to_milliseconds(
            self.update_timestamp)
        if record.id is not None:
            variant.names.extend(record.id.split(';'))
        variant.start = record.start          # 0-based inclusive
        variant.end = record.stop             # 0-based exclusive
        variant.reference_bases = record.ref
        if record.alts is not None:
            variant.alternate_bases.extend(list(record.alts))
        # record.filter and record.qual are also available, when supported
        # by GAVariant.
        for key, value in record.info.iteritems():
            if value is not None:
                if isinstance(value, str):
                    value = value.split(',')
                variant.info[key].values.extend(_encodeValue(value))
        for call_set in call_sets:
            pysam_call = record.samples[call_set.sample_id.encode()]
            ga_call = variant.calls.add()
            ga_call.call_set_id = str(call_set.id)
            ga_call.call_set_name = str(call_set.sample_id)
            self._update_ga_call(ga_call, pysam_call)
        variant.id = self.get_variant_id(variant)
        return variant

    # TODO refactor the common code in these external methods.
    def run_search(self, request, call_sets, response_builder):
        # First get the VCF file for this reference.
        vcf_file = self.vcf_files.filter(
            VcfFile.reference_name == request.reference_name).first()
        if vcf_file is not None:
            var_file = self.get_file_handle(
                pysam.VariantFile, vcf_file.data_url,
                index_filename=vcf_file.index_url)
            the_start = request.start
            if request.page_token:
                compound_id = registry.VariantCompoundId.parse(
                    request.page_token)
                the_start = int(compound_id.start)
            reference_name, start, end = self.sanitizeVariantFileFetch(
                request.reference_name, the_start, request.end)
            cursor = var_file.fetch(reference_name, start, end)
            variant = None
            if request.page_token:
                # If we have a pagetoken, we must skip ahead until we reach the
                # variant with this ID.
                for record in cursor:
                    variant = self.convert_variant(record, call_sets)
                    if variant.id == request.page_token:
                        break
            for record in cursor:
                variant = self.convert_variant(record, call_sets)
                response_builder.addValue(variant)
                if response_builder.isFull():
                    break
            if response_builder.isFull() and not is_empty_iter(cursor):
                response_builder.setNextPageToken(variant.id)

    def run_get_variant(self, compound_id):
        vcf_file = self.vcf_files.filter(
            VcfFile.reference_name == compound_id.reference_name).first()
        if vcf_file is None:
            raise exceptions.ObjectNotFoundException(compound_id)
        var_file = self.get_file_handle(
            pysam.VariantFile, vcf_file.data_url,
            index_filename=vcf_file.index_url)
        start = int(compound_id.start)
        reference_name, start, end = self.sanitizeVariantFileFetch(
                compound_id.reference_name, start, start + 1)
        cursor = var_file.fetch(reference_name, start, end)
        id_str = str(compound_id)
        for record in cursor:
            variant = self.convert_variant(record, self.call_sets)
            if variant.id == id_str:
                return variant
            elif record.start > start:
                raise exceptions.ObjectNotFoundException()
        raise exceptions.ObjectNotFoundException(compound_id)

    def get_variants(self, reference_name, start, end, call_set_ids=[]):
        """
        Convenience method for testing. Returns an iterator over the
        specified set of variants.

        TODO remove this method once the corresponding data driven testing
        code has been refactored to use the main server pathway and local
        client interface.
        """
        # Get the list of call sets. This algorithm is inefficient, but it
        # doesn't matter as this method is only used in testing.
        call_sets = []
        for call_set in self.call_sets:
            if str(call_set.id) in call_set_ids:
                call_sets.append(call_set)
        vcf_file = self.vcf_files.filter(
            VcfFile.reference_name == reference_name).first()
        if vcf_file is not None:
            var_file = self.get_file_handle(
                pysam.VariantFile, vcf_file.data_url,
                index_filename=vcf_file.index_url)
            cursor = var_file.fetch(reference_name.encode(), start, end)
            for record in cursor:
                yield self.convert_variant(record, call_sets)


##################################################
#
# Reads
#
##################################################


class HtslibReadGroupSet(registry.ReadGroupSet, PysamMixin):

    __tablename__ = 'HtslibReadGroupSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReadGroupSet.id'),
        primary_key=True)
    data_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    index_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    assembly_identifier = sqlalchemy.Column(sqlalchemy.String)
    """
    The Assembly Identifier read from the BAM file header.
    """

    __mapper_args__ = {
        'polymorphic_identity': 'HtslibReadGroupSet',
    }

    def __init__(self, name, data_url, index_url):
        super(HtslibReadGroupSet, self).__init__(name)
        self.data_url = data_url
        self.index_url = index_url
        alignment_file = self.__open_file(data_url, index_url)
        try:
            self.__read_alignment_file(alignment_file)
        finally:
            alignment_file.close()

    def check_references(self, reference_set):
        """
        Ensure that references are consistent between the DB reference set
        and the references found in the BAM header.
        """
        alignment_file = self.__open_file(self.data_url, self.index_url)
        try:
            # Get the reference names from the alignment file.
            header_reference_names = set()
            for reference_name in alignment_file.references:
                header_reference_names.add(reference_name)
            db_reference_names = set(
                reference.name.encode() for reference
                in reference_set.references)
            if not db_reference_names.issubset(header_reference_names):
                missing = db_reference_names - header_reference_names
                raise exceptions.ReadGroupReferenceNotFound(
                    self.data_url, sorted(missing), reference_set.name)
        finally:
            alignment_file.close()

    def __get_read_group_id(self, name):
        """
        Returns the read group ID associated with the specified name. The
        name must be bytes encoded _not_ unicode. If the name is not known
        we return the empty string.
        """
        # TODO is returning the empty string the correct behaviour here? The
        # read doesn't belong to a read group in the file, so should we have a
        # synthetic read group to catch those that aren't explicitlu assigned
        # to one??
        # TODO avoid creating the dictionary here.
        name_map = {rg.name.encode(): str(rg.id) for rg in self.read_groups}
        return name_map.get(name, "")

    def __open_file(self, data_url, index_url):
        # We need to check to see if the path exists here as pysam does
        # not throw an error if the index is missing.
        if not os.path.exists(index_url):
            raise exceptions.FileOpenFailedException(index_url)
        try:
            return pysam.AlignmentFile(data_url, filepath_index=index_url)
        except IOError as exception:
            # IOError thrown when the index file passed in is not actually
            # an index file... may also happen in other cases?
            raise exceptions.DataException(exception.message)

    def __create_read_group(self, rg_header):
        """
        Returns a new read group for the specified SAM RG rg_header.
        """
        read_group = registry.ReadGroup(rg_header['ID'])
        read_group.sample_id = rg_header.get('SM', "")
        read_group.description = rg_header.get('DS', "")
        read_group.predicted_insert_size = int(rg_header.get('PI', 0))
        read_group.experiment = registry.Experiment(
            instrument_model=rg_header.get('PL', ""),
            sequencing_center=rg_header.get('CN', ""),
            description=rg_header.get('DS', ""),
            library=rg_header.get('LB', ""),
            platform_unit=rg_header.get('PU', ""),
            run_time=rg_header.get('DT', "")
        )
        return read_group

    def __read_alignment_file(self, alignment_file):
        """
        Populates the instance variables of this read group set
        based on the specified pysam AlignmentFile instance.
        """
        header = alignment_file.header
        # Create the read groups.
        if 'RG' not in header or len(header['RG']) == 0:
            # If there is no read group header information, we assume there
            # is a single read group in the file, and give the default name
            read_group = registry.ReadGroup(
                name=DEFAULT_READGROUP_NAME, sample_id="")
            self.read_groups.append(read_group)
        else:
            for rg_header in header['RG']:
                self.read_groups.append(self.__create_read_group(rg_header))
        # Create the Programs for this ReadGroupSet.
        if 'PG' in header:
            for htslib_program in header['PG']:
                program = registry.Program()
                program.command_line = htslib_program.get('CL', "")
                program.name = htslib_program.get('PN', "")
                # TODO when we can get the PP values from the header, make a
                # map # from the SAM IDs to our external IDs and map this
                # to previous program ID.
                # htslib_id = htslibProgram['ID']
                # htslib_prev_if  = htslibProgram.get('PP', "")
                program.version = htslib_program.get('VN', "")
                self.programs.append(program)
        # Get the read stats. We cannot currently get the total number
        # of bases, as pysam does not provide it.
        self.read_stats = registry.ReadStats(
            aligned_read_count=alignment_file.mapped,
            unaligned_read_count=alignment_file.unmapped)
        # Get the reference information from the header.
        self.assembly_identifier = None
        for sq_header in header['SQ']:
            if 'AS' not in sq_header:
                sq_header = parseMalformedBamHeader(sq_header)
            name = sq_header.get('AS', None)
            if self.assembly_identifier is None:
                self.assembly_identifier = name
            elif self.assembly_identifier != name:
                raise exceptions.MultipleReferenceSetsInReadGroupSet(
                    self.data_url, name, self.assembly_identifier)

    def convertReadAlignment(self, samFile, read, reference):
        """
        Convert a pysam ReadAlignment to a GA4GH ReadAlignment
        """
        # TODO fill out remaining fields
        # TODO refine in tandem with code in converters module
        ret = protocol.ReadAlignment()
        # ret.fragmentId = 'TODO'
        ret.aligned_quality.extend(read.query_qualities)
        ret.aligned_sequence = read.query_sequence
        if SamFlags.isFlagSet(read.flag, SamFlags.READ_UNMAPPED):
            ret.ClearField("alignment")
        else:
            ret.alignment.CopyFrom(protocol.LinearAlignment())
            ret.alignment.mapping_quality = read.mapping_quality
            ret.alignment.position.CopyFrom(protocol.Position())
            ret.alignment.position.reference_name = samFile.getrname(
                read.reference_id)
            ret.alignment.position.position = read.reference_start
            ret.alignment.position.strand = protocol.POS_STRAND
            if SamFlags.isFlagSet(read.flag, SamFlags.READ_REVERSE_STRAND):
                ret.alignment.position.strand = protocol.NEG_STRAND
            for operation, length in read.cigar:
                gaCigarUnit = ret.alignment.cigar.add()
                gaCigarUnit.operation = SamCigar.int2ga(operation)
                gaCigarUnit.operation_length = length
                gaCigarUnit.reference_sequence = ""  # TODO fix this!
        ret.duplicate_fragment = SamFlags.isFlagSet(
            read.flag, SamFlags.DUPLICATE_READ)
        ret.failed_vendor_quality_checks = SamFlags.isFlagSet(
            read.flag, SamFlags.FAILED_QUALITY_CHECK)
        ret.fragment_length = read.template_length
        ret.fragment_name = read.query_name
        # TODO default read group name
        read_group_name = DEFAULT_READGROUP_NAME
        for key, value in read.tags:
            if key == b"RG":
                read_group_name = value
            else:
                ret.info[key].values.add().string_value = str(value)
        ret.read_group_id = self.__get_read_group_id(read_group_name)
        if SamFlags.isFlagSet(read.flag, SamFlags.MATE_UNMAPPED):
            ret.next_mate_position.Clear()
        else:
            ret.next_mate_position.Clear()
            if read.next_reference_id != -1:
                ret.next_mate_position.reference_name = samFile.getrname(
                    read.next_reference_id)
            else:
                ret.next_mate_position.reference_name = ""
            ret.next_mate_position.position = read.next_reference_start
            ret.next_mate_position.strand = protocol.POS_STRAND
            if SamFlags.isFlagSet(read.flag, SamFlags.MATE_REVERSE_STRAND):
                ret.next_mate_position.strand = protocol.NEG_STRAND
        if SamFlags.isFlagSet(read.flag, SamFlags.READ_PAIRED):
            ret.number_reads = 2
        else:
            ret.number_reads = 1
        ret.read_number = -1
        if SamFlags.isFlagSet(read.flag, SamFlags.FIRST_IN_PAIR):
            if SamFlags.isFlagSet(read.flag, SamFlags.SECOND_IN_PAIR):
                ret.read_number = 2
            else:
                ret.read_number = 0
        elif SamFlags.isFlagSet(read.flag, SamFlags.SECOND_IN_PAIR):
            ret.read_number = 1
        ret.improper_placement = not SamFlags.isFlagSet(
            read.flag, SamFlags.READ_PROPER_PAIR)
        ret.secondary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SECONDARY_ALIGNMENT)
        ret.supplementary_alignment = SamFlags.isFlagSet(
            read.flag, SamFlags.SUPPLEMENTARY_ALIGNMENT)
        # TODO is this really a unique ID for this particular read? Should we
        # not hash something as well to be certain?
        compound_id = registry.ReadAlignmentCompoundId(
            None, ret.read_group_id, str(reference.id),
            str(ret.alignment.position.position), ret.fragment_name)
        ret.id = str(compound_id)
        return ret

    def run_search(self, request, reference, read_groups, response_builder):
        start = request.start
        end = request.end
        reference_name = reference.name.encode()
        if request.page_token:
            compound_id = registry.ReadAlignmentCompoundId.parse(
                request.page_token)
            start = int(compound_id.position)
        start, end = self.sanitizeAlignmentFileFetch(start, end)
        alignment_file = self.__open_file(self.data_url, self.index_url)
        read_alignments = alignment_file.fetch(reference_name, start, end)
        if request.page_token:
            # Skip ahead until we see the read with this ID.
            for alignment in read_alignments:
                ga_alignment = self.convertReadAlignment(
                    alignment_file, alignment, reference)
                if ga_alignment.id == request.page_token:
                    break
        for alignment in read_alignments:
            ga_alignment = self.convertReadAlignment(
                alignment_file, alignment, reference)
            if ga_alignment.read_group_id in request.read_group_ids:
                response_builder.addValue(ga_alignment)
            if response_builder.isFull():
                break
        if response_builder.isFull() and not is_empty_iter(read_alignments):
            response_builder.setNextPageToken(ga_alignment.id)


class PysamFileHandleCache(object):
    """
    Cache for opened file handles. We use a deque which has the
    advantage to have push/pop operations in O(1) We always add
    elements on the left of the deque and pop elements from the right.
    When a file is accessed via getFileHandle, its priority gets
    updated, it is put at the "top" of the deque.
    """

    def __init__(self):
        self._cache = collections.deque()
        self._memoTable = dict()
        # Initialize the value even if it will be set up by the config
        self._maxCacheSize = 50

    def setMaxCacheSize(self, size):
        """
        Sets the maximum size of the cache
        """
        if size <= 0:
            raise ValueError(
                "The size of the cache must be a strictly positive value")
        self._maxCacheSize = size

    def _add(self, dataFile, handle):
        """
        Add a file handle to the left of the deque
        """
        self._cache.appendleft((dataFile, handle))

    def _update(self, dataFile, handle):
        """
        Update the priority of the file handle. The element is first
        removed and then added to the left of the deque.
        """
        self._cache.remove((dataFile, handle))
        self._add(dataFile, handle)

    def _removeLru(self):
        """
        Remove the least recently used file handle from the cache.
        The pop method removes an element from the right of the deque.
        Returns the name of the file that has been removed.
        """
        (dataFile, handle) = self._cache.pop()
        handle.close()
        return dataFile

    def getCachedFiles(self):
        """
        Returns all file names stored in the cache.
        """
        return self._memoTable.keys()

    def get_file_handle(self, url, open_func, *args, **kwargs):
        """
        Returns handle associated to the filename. If the file is
        already opened, update its priority in the cache and return
        its handle. Otherwise, open the file using openMethod, store
        it in the cache and return the corresponding handle.
        """
        if url in self._memoTable:
            handle = self._memoTable[url]
            self._update(url, handle)
            return handle
        else:
            try:
                handle = open_func(*args, **kwargs)
            except ValueError:
                raise exceptions.FileOpenFailedException(url)

            self._memoTable[url] = handle
            self._add(url, handle)
            if len(self._memoTable) > self._maxCacheSize:
                url = self._removeLru()
                del self._memoTable[url]
            return handle


# LRU cache of open file handles
file_handle_cache = PysamFileHandleCache()
