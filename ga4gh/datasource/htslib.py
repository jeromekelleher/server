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

import ga4gh.registry as registry

_nothing = object()


def is_empty_iter(it):
    """
    Return True iff the iterator is empty or exhausted
    """
    return next(it, _nothing) is _nothing



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


class HtslibVariantSet(registry.VariantSet):
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

        # self._updateMetadata(varFile)
        # self._updateCallSetIds(varFile)
        # self._updateVariantAnnotationSets(varFile, dataUrl)


    def run_search(self, request, response_builder):
        print("running search on ", request, response_builder)
        print("reference name =", type(request.reference_name), request.reference_name)
        # First get the VCF file for this reference.
        vcf_file = self.vcf_files.filter(
            VcfFile.reference_name == request.reference_name).first()
        if vcf_file is not None:
            print("vcf_file", vcf_file.data_url)
        else:
            print("Not found")


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

    def getFileHandle(self, dataFile, openMethod):
        """
        Returns handle associated to the filename. If the file is
        already opened, update its priority in the cache and return
        its handle. Otherwise, open the file using openMethod, store
        it in the cache and return the corresponding handle.
        """
        if dataFile in self._memoTable:
            handle = self._memoTable[dataFile]
            self._update(dataFile, handle)
            return handle
        else:
            try:
                handle = openMethod(dataFile)
            except ValueError:
                raise exceptions.FileOpenFailedException(dataFile)

            self._memoTable[dataFile] = handle
            self._add(dataFile, handle)
            if len(self._memoTable) > self._maxCacheSize:
                dataFile = self._removeLru()
                del self._memoTable[dataFile]
            return handle


# LRU cache of open file handles
fileHandleCache = PysamFileHandleCache()



