"""
The GA4GH data model. Defines all the methods required to translate
data in existing formats into GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import glob
import os
import json
import tempfile
import shutil
import atexit

import ga4gh.exceptions as exceptions


def _cleanupHtslibsMess(indexDir):
    """
    Cleanup the mess that htslib has left behind with the index files.
    This is a temporary measure until we get a good interface for
    dealing with indexes for remote files.
    """
    if os.path.exists(indexDir):
        shutil.rmtree(indexDir)


class DatamodelObject(object):
    """
    Superclass of all datamodel types
    """
    def __init__(self):
        # TODO move common functionality into this class from subclasses
        pass


class PysamDatamodelMixin(object):
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
    def sanitizeAlignmentFileFetch(
            cls, reference_name=None, start=None, end=None):
        if reference_name is not None:
            reference_name = cls.sanitizeString(
                reference_name, 'reference_name')
        if start is not None:
            start = cls.sanitizeInt(
                start, cls.samMin, cls.samMaxStart, 'start')
        if end is not None:
            end = cls.sanitizeInt(end, cls.samMin, cls.samMaxEnd, 'end')
        if start is not None and end is not None:
            cls.assertValidRange(start, end, 'start', 'end')
        return reference_name, start, end

    @classmethod
    def sanitizeFastaFileFetch(cls, start=None, end=None):
        if start is not None:
            start = cls.sanitizeInt(
                start, cls.fastaMin, cls.fastaMax, 'start')
        if end is not None:
            end = cls.sanitizeInt(
                end, cls.fastaMin, cls.fastaMax, 'end')
        if start is not None and end is not None:
            cls.assertValidRange(start, end, 'start', 'end')
        return start, end

    @classmethod
    def sanitizeGetRName(cls, reference_id):
        cls.assertInt(reference_id, 'reference_id')
        cls.assertInRange(
            reference_id, cls.rNameMin, cls.rNameMax, 'reference_id')

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
        if not (isinstance(attr, int) or isinstance(attr, long)):
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

    def _setAccessTimes(self, directoryPath):
        """
        Sets the creationTime and accessTime for this file system based
        DatamodelObject. This is derived from the ctime of the specified
        directoryPath.
        """
        ctimeInMillis = int(os.path.getctime(directoryPath) * 1000)
        self._creationTime = ctimeInMillis
        self._updatedTime = ctimeInMillis

    def _scanDataFiles(self, dataDir, patterns):
        """
        Scans the specified directory for files with the specified globbing
        pattern and calls self._addDataFile for each. Raises an
        EmptyDirException if no data files are found.
        """
        numDataFiles = 0
        for pattern in patterns:
            scanPath = os.path.join(dataDir, pattern)
            for filename in glob.glob(scanPath):
                self._addDataFile(filename)
                numDataFiles += 1
        # This is a temporary workaround to allow us to use htslib's
        # facility for working with remote files. The urls.json is
        # definitely not a good idea and will be replaced later.
        # We make a temporary file for each process so that it
        # downloads its own copy and we are sure it's not overwriting
        # the copy of another process. We then register a cleanup
        # handler to get rid of these files on exit.
        urlSource = os.path.join(dataDir, "urls.json")
        if os.path.exists(urlSource):
            with open(urlSource) as jsonFile:
                urls = json.load(jsonFile)["urls"]
            indexDir = tempfile.mkdtemp(prefix="htslib_mess.")
            cwd = os.getcwd()
            os.chdir(indexDir)
            for url in urls:
                self._addDataFile(url)
                numDataFiles += 1
            os.chdir(cwd)
            atexit.register(_cleanupHtslibsMess, indexDir)
        if numDataFiles == 0:
            raise exceptions.EmptyDirException(dataDir, patterns)
