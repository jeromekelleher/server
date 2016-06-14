"""
Module defining concrete instances of the top level containers using
standard data files interpreted using pysam/htslib.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import hashlib
import collections

import pysam
import sqlalchemy
import sqlalchemy.orm as orm
import google.protobuf.struct_pb2 as struct_pb2

import ga4gh.registry as registry
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions

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


class HtslibReferenceSet(registry.ReferenceSet):
    __tablename__ = 'HtslibReferenceSet'
    id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey('ReferenceSet.id'),
        primary_key=True)
    fasta_url = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    __mapper_args__ = {
        'polymorphic_identity':'HtslibReferenceSet',
    }

    def __init__(self, name, fasta_url):
        super(HtslibReferenceSet, self).__init__(name)
        self.fasta_url = fasta_url
        with pysam.FastaFile(fasta_url) as fasta_file:
            for reference_name in fasta_file.references:
                reference = registry.Reference(reference_name)
                bases = fasta_file.fetch(reference_name)
                reference.md5checksum = hashlib.md5(bases).hexdigest()
                reference.length = len(bases)
                self.references.append(reference)
        # Now that we have all the references, set the md5checksum
        self._compute_md5checksum()


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
        'polymorphic_identity':'HtslibVariantSet',
    }

    def __init__(self, name, data_urls, index_urls):
        super(HtslibVariantSet, self).__init__(name)
        print("creating htslib variant set for ", data_urls, index_urls)
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
            if reference_name in self.vcf_files:
                raise exceptions.OverlappingVcfException(
                    vcf_file.data_url, reference_name)


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
                # VCF because we need to find the file that maps to a particular
                # reference.
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
            var_file = file_handle_cache.get_file_handle(
                vcf_file.data_url, pysam.VariantFile, vcf_file.data_url,
                index_filename=vcf_file.index_url)
            the_start = request.start
            skip_required = False
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
            if response_builder.isFull():
                response_builder.setNextPageToken(variant.id)

    def run_get_variant(self, compound_id):
        vcf_file = self.vcf_files.filter(
            VcfFile.reference_name == compound_id.reference_name).first()
        if vcf_file is None:
            raise exceptions.ObjectNotFoundException(compound_id)
        var_file = file_handle_cache.get_file_handle(
            vcf_file.data_url, pysam.VariantFile, vcf_file.data_url,
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
        raise exceptions.ObjectNotFoundException(compoundId)


    def getVcfHeaderReferenceSetName(self):
        """
        Returns the name of the reference set from the VCF header.
        """
        # TODO implemenent
        return None


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


