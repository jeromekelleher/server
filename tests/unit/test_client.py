"""
Tests for the client
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import mock

import ga4gh.protocol as protocol
import ga4gh.backend as backend
import ga4gh.client as client
import tests.utils as utils
import ga4gh.exceptions as exceptions


class TestSearchMethodsCallRunRequest(unittest.TestCase):
    """
    Test that search methods call lower-level functionality correctly
    """
    def setUp(self):
        self.httpClient = client.HttpClient("http://example.com")
        self.httpClient._runSearchRequest = mock.Mock()
        self.httpClient._runGetRequest = mock.Mock()
        self.objectId = "SomeId"
        self.objectName = "objectName"
        self.datasetId = "datasetId"
        self.variantSetId = "variantSetId"
        self.variantAnnotationSetId = "variantAnnotationSetId"
        self.featureSetId = "featureSetId"
        self.parentId = "parentId"
        self.feature = "feature"
        self.referenceSetId = "referenceSetId"
        self.referenceId = "referenceId"
        self.readGroupIds = ["readGroupId"]
        self.referenceName = "referenceName"
        self.bioSampleId = "bioSampleId"
        self.bioSampleName = "bioSampleName"
        self.individualName = "individualName"
        self.individualId = "individualId"
        self.start = 100
        self.end = 101
        self.referenceName = "referenceName"
        self.callSetIds = ["id1", "id2"]
        self.pageSize = 1000
        self.httpClient.setPageSize(self.pageSize)
        self.assemblyId = "assemblyId"
        self.accession = "accession"
        self.md5checksum = "md5checksum"

    def testSetPageSize(self):
        testClient = client.AbstractClient()
        # pageSize is None by default
        self.assertIsNone(testClient.getPageSize())
        for pageSize in [1, 10, 100]:
            testClient.setPageSize(pageSize)
            self.assertEqual(testClient.getPageSize(), pageSize)

    def testSearchVariants(self):
        request = protocol.SearchVariantsRequest()
        request.reference_name = self.referenceName
        request.start = self.start
        request.end = self.end
        request.variant_set_id = self.variantSetId
        request.call_set_ids.extend(self.callSetIds)
        request.page_size = self.pageSize
        self.httpClient.searchVariants(
            self.variantSetId, start=self.start, end=self.end,
            referenceName=self.referenceName, callSetIds=self.callSetIds)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "variants", protocol.SearchVariantsResponse)

    def testSearchDatasets(self):
        request = protocol.SearchDatasetsRequest()
        request.page_size = self.pageSize
        self.httpClient.searchDatasets()
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "datasets", protocol.SearchDatasetsResponse)

    def testSearchVariantSets(self):
        request = protocol.SearchVariantSetsRequest()
        request.dataset_id = self.datasetId
        request.page_size = self.pageSize
        self.httpClient.searchVariantSets(self.datasetId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "variantsets", protocol.SearchVariantSetsResponse)

    def testSearchVariantAnnotationSets(self):
        request = protocol.SearchVariantAnnotationSetsRequest()
        request.variant_set_id = self.variantSetId
        request.page_size = self.pageSize
        self.httpClient.searchVariantAnnotationSets(self.variantSetId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "variantannotationsets",
            protocol.SearchVariantAnnotationSetsResponse)

    def testSearchVariantAnnotations(self):
        request = protocol.SearchVariantAnnotationsRequest()
        request.variant_annotation_set_id = self.variantAnnotationSetId
        request.page_size = self.pageSize
        request.reference_name = self.referenceName
        request.reference_id = self.referenceId
        request.start = self.start
        request.end = self.end
        self.httpClient.searchVariantAnnotations(
            self.variantAnnotationSetId,
            referenceName=self.referenceName,
            start=self.start,
            end=self.end,
            effects=[],
            referenceId=self.referenceId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "variantannotations",
            protocol.SearchVariantAnnotationsResponse)
        with self.assertRaises(exceptions.BadRequestException):
            self.httpClient.searchVariantAnnotations(
                self.variantAnnotationSetId,
                referenceName=self.referenceName,
                start=self.start,
                end=self.end,
                effects=[{"term": "just a term"}, {"id": "an id"}],
                referenceId=self.referenceId)

    def testSearchFeatures(self):
        request = protocol.SearchFeaturesRequest()
        request.feature_set_id = self.featureSetId
        request.parent_id = self.parentId
        request.page_size = self.pageSize
        request.reference_name = self.referenceName
        request.start = self.start
        request.end = self.end
        request.feature_types.append(self.feature)
        self.httpClient.searchFeatures(
            self.featureSetId, parentId=self.parentId,
            referenceName=self.referenceName, start=self.start,
            end=self.end, featureTypes=[self.feature])
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "features", protocol.SearchFeaturesResponse)

    def testSearchFeatureSets(self):
        request = protocol.SearchFeatureSetsRequest()
        request.dataset_id = self.datasetId
        request.page_size = self.pageSize
        self.httpClient.searchFeatureSets(self.datasetId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "featuresets", protocol.SearchFeatureSetsResponse)

    def testSearchReferenceSets(self):
        request = protocol.SearchReferenceSetsRequest()
        request.page_size = self.pageSize
        request.accession = self.accession
        request.md5checksum = self.md5checksum
        request.assembly_id = self.assemblyId
        self.httpClient.searchReferenceSets(
            accession=self.accession, md5checksum=self.md5checksum,
            assemblyId=self.assemblyId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "referencesets", protocol.SearchReferenceSetsResponse)

    def testSearchReferences(self):
        request = protocol.SearchReferencesRequest()
        request.reference_set_id = self.referenceSetId
        request.page_size = self.pageSize
        request.accession = self.accession
        request.md5checksum = self.md5checksum
        self.httpClient.searchReferences(
            self.referenceSetId, accession=self.accession,
            md5checksum=self.md5checksum)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "references", protocol.SearchReferencesResponse)

    def testSearchReadGroupSets(self):
        request = protocol.SearchReadGroupSetsRequest()
        request.dataset_id = self.datasetId
        request.name = self.objectName
        request.bio_sample_id = self.bioSampleId
        request.page_size = self.pageSize
        self.httpClient.searchReadGroupSets(
            self.datasetId,
            name=self.objectName,
            bioSampleId=self.bioSampleId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "readgroupsets", protocol.SearchReadGroupSetsResponse)

    def testSearchCallSets(self):
        request = protocol.SearchCallSetsRequest()
        request.variant_set_id = self.variantSetId
        request.name = self.objectName
        request.bio_sample_id = self.bioSampleId
        request.page_size = self.pageSize
        self.httpClient.searchCallSets(
            self.variantSetId,
            name=self.objectName,
            bioSampleId=self.bioSampleId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "callsets", protocol.SearchCallSetsResponse)

    def testSearchReads(self):
        request = protocol.SearchReadsRequest()
        request.read_group_ids.extend(self.readGroupIds)
        request.reference_id = self.referenceId
        request.start = self.start
        request.end = self.end
        request.page_size = self.pageSize
        self.httpClient.searchReads(
            self.readGroupIds, referenceId=self.referenceId,
            start=self.start, end=self.end)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "reads", protocol.SearchReadsResponse)

    def testSearchBioSamples(self):
        request = protocol.SearchBioSamplesRequest()
        request.dataset_id = self.datasetId
        request.name = self.bioSampleName
        request.individual_id = self.individualId
        request.page_size = self.pageSize
        self.httpClient.searchBioSamples(
            self.datasetId, self.bioSampleName, self.individualId)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "biosamples", protocol.SearchBioSamplesResponse)

    def testSearchIndividuals(self):
        request = protocol.SearchIndividualsRequest()
        request.dataset_id = self.datasetId
        request.name = self.individualName
        request.page_size = self.pageSize
        self.httpClient.searchIndividuals(
            self.datasetId, self.individualName)
        self.httpClient._runSearchRequest.assert_called_once_with(
            request, "individuals", protocol.SearchIndividualsResponse)

    def testGetReferenceSet(self):
        self.httpClient.getReferenceSet(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "referencesets", protocol.ReferenceSet, self.objectId)

    def testGetVariantAnnotationSet(self):
        self.httpClient.getVariantAnnotationSet(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "variantannotationsets", protocol.VariantAnnotationSet,
            self.objectId)

    def testGetVariantSet(self):
        self.httpClient.getVariantSet(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "variantsets", protocol.VariantSet, self.objectId)

    def testGetReference(self):
        self.httpClient.getReference(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "references", protocol.Reference, self.objectId)

    def testGetReadGroupSets(self):
        self.httpClient.getReadGroupSet(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "readgroupsets", protocol.ReadGroupSet, self.objectId)

    def testGetReadGroup(self):
        self.httpClient.getReadGroup(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "readgroups", protocol.ReadGroup, self.objectId)

    def testGetCallSets(self):
        self.httpClient.getCallSet(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "callsets", protocol.CallSet, self.objectId)

    def testGetDatasets(self):
        self.httpClient.getDataset(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "datasets", protocol.Dataset, self.objectId)

    def testGetVariant(self):
        self.httpClient.getVariant(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "variants", protocol.Variant, self.objectId)

    def testGetBioSample(self):
        self.httpClient.getBioSample(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "biosamples", protocol.BioSample, self.objectId)

    def testGetIndividual(self):
        self.httpClient.getIndividual(self.objectId)
        self.httpClient._runGetRequest.assert_called_once_with(
            "individuals", protocol.Individual, self.objectId)


class DatamodelObjectWrapper(object):
    """
    Thin wrapper class that allows us to treat data model objects uniformly.
    We should update the data model interface so that objects are always
    returned so that we always call get_protobuf on the results.
    Variants and Reads are the exceptions here.
    """
    def __init__(self, gaObject):
        self.gaObject = gaObject

    def get_protobuf(self):
        return self.gaObject


class DummyResponse(object):
    """
    Stand in for requests Response object;
    """
    def __init__(self, text):
        self.text = text
        self.status_code = 200


class DummyRequestsSession(object):
    """
    Takes the place of a requests session so that we can check that all
    values are sent and received correctly.
    """
    def __init__(self, backend, urlPrefix):
        self._backend = backend
        self._urlPrefix = urlPrefix
        self._getMethodMap = {
            "datasets": self._backend.runGetDataset,
            "referencesets": self._backend.runGetReferenceSet,
            "references": self._backend.runGetReference,
            "variantsets": self._backend.runGetVariantSet,
            "variants": self._backend.runGetVariant,
            "readgroupsets": self._backend.runGetReadGroupSet,
            "readgroups": self._backend.runGetReadGroup,
        }
        self._searchMethodMap = {
            "datasets": self._backend.runSearchDatasets,
            "referencesets": self._backend.runSearchReferenceSets,
            "references": self._backend.runSearchReferences,
            "variantsets": self._backend.runSearchVariantSets,
            "variants": self._backend.runSearchVariants,
            "readgroupsets": self._backend.runSearchReadGroupSets,
            "reads": self._backend.runSearchReads,
        }
        self.headers = {}

    def checkSessionParameters(self):
        contentType = "Content-type"
        assert contentType in self.headers
        assert self.headers[contentType] == "application/json"

    def get(self, url, params):
        # TODO add some more checks for params to see if Key is set,
        # and we're not sending any extra stuff.
        self.checkSessionParameters()
        assert url.startswith(self._urlPrefix)
        suffix = url[len(self._urlPrefix):]
        basesSuffix = "/bases"
        splits = suffix.split("/")
        if suffix.endswith(basesSuffix):
            # ListReferenceBases is an oddball and needs to be treated
            # separately.
            assert splits[0] == ''
            assert splits[1] == 'references'
            id_ = splits[2]
            assert splits[3] == 'bases'
            # This is all very ugly --- see the comments in the LocalClient
            # for why we need to do this. Definitely needs to be fixed.
            args = dict(params)
            if args[u'end'] == u'0':
                del args['end']
            if args['pageToken'] is "":
                del args['pageToken']
            result = self._backend.runListReferenceBases(id_, args)
        else:
            assert len(splits) == 3
            assert splits[0] == ''
            datatype, id_ = splits[1:]
            assert datatype in self._getMethodMap
            method = self._getMethodMap[datatype]
            result = method(id_)
        return DummyResponse(result)

    def post(self, url, params=None, data=None):
        self.checkSessionParameters()
        assert url.startswith(self._urlPrefix)
        suffix = url[len(self._urlPrefix):]
        searchSuffix = "/search"
        assert suffix.startswith("/")
        assert suffix.endswith(searchSuffix)
        datatype = suffix[1:-len(searchSuffix)]
        assert datatype in self._searchMethodMap
        method = self._searchMethodMap[datatype]
        result = method(data)
        return DummyResponse(result)


class DummyHttpClient(client.HttpClient):
    """
    Client in which we intercept calls to the underlying requests connection.
    """
    def __init__(self, backend):
        self._urlPrefix = "http://example.com"
        super(DummyHttpClient, self).__init__(self._urlPrefix)
        self._session = DummyRequestsSession(backend, self._urlPrefix)
        self._setupHttpSession()


class ExhaustiveListingsMixin(object):
    """
    Tests exhaustive listings using the high-level API with a Simulated
    backend.
    """
    @classmethod
    def setUpClass(cls):
        registry_db = utils.create_simulated_registry_db(
            random_seed=100, num_datasets=3,
            num_variant_sets=3, num_calls=3,
            num_reference_sets=3, num_references_per_reference_set=3,
            num_read_group_sets=3, num_read_groups_per_read_group_set=3)
        cls.backend = backend.Backend(registry_db)
        cls.registry_db = registry_db

    def setUp(self):
        self.client = self.getClient()

    def verifyObjectList(self, gaObjects, datamodelObjects, getMethod):
        """
        Verifies that the specified list of protocol objects corresponds
        to the specified list of datamodel objects.
        """
        num_objects = 0
        for gaObject, datamodelObject in utils.zipLists(
                gaObjects, datamodelObjects):
            self.assertEqual(gaObject, datamodelObject.get_protobuf())
            otherGaObject = getMethod(gaObject.id)
            self.assertEqual(gaObject, otherGaObject)
            num_objects += 1
        self.assertGreater(num_objects, 0)

    def testAllDatasets(self):
        datasets = list(self.client.searchDatasets())
        self.verifyObjectList(
            datasets, self.registry_db.get_datasets(), self.client.getDataset)

    def testAllReferenceSets(self):
        referenceSets = list(self.client.searchReferenceSets())
        self.verifyObjectList(
            referenceSets, self.registry_db.get_reference_sets(),
            self.client.getReferenceSet)

    def testAllReferences(self):
        for referenceSet in self.client.searchReferenceSets():
            references = list(self.client.searchReferences(referenceSet.id))
            datamodelReferences = self.registry_db.get_reference_set(
                referenceSet.id).references
            self.verifyObjectList(
                references, datamodelReferences, self.client.getReference)
            for datamodelReference in datamodelReferences:
                bases = self.client.listReferenceBases(datamodelReference.id)
                otherBases = datamodelReference.run_get_bases(
                    0, datamodelReference.length)
                self.assertEqual(bases, otherBases)

    def testAllVariantSets(self):
        for dataset in self.client.searchDatasets():
            variantSets = list(self.client.searchVariantSets(dataset.id))
            datamodelVariantSets = list(self.registry_db.get_dataset(
                dataset.id).variant_sets)
            self.verifyObjectList(
                variantSets, datamodelVariantSets, self.client.getVariantSet)

    def testAllReadGroupSets(self):
        for dataset in self.client.searchDatasets():
            readGroupSets = list(self.client.searchReadGroupSets(dataset.id))
            datamodelReadGroupSets = list(self.registry_db.get_dataset(
                dataset.id).read_group_sets)
            self.verifyObjectList(
                readGroupSets, datamodelReadGroupSets,
                self.client.getReadGroupSet)
            # Check the readGroups.
            for readGroupSet, datamodelReadGroupSet in zip(
                    readGroupSets, datamodelReadGroupSets):
                datamodelReadGroups = datamodelReadGroupSet.read_groups
                self.verifyObjectList(
                    readGroupSet.read_groups, datamodelReadGroups,
                    self.client.getReadGroup)


class TestExhaustiveListingsHttp(ExhaustiveListingsMixin, unittest.TestCase):
    """
    Tests the exhaustive listings using the HTTP client.
    """

    def getClient(self):
        return DummyHttpClient(self.backend)


class TestExhaustiveListingsLocal(ExhaustiveListingsMixin, unittest.TestCase):
    """
    Tests the exhaustive listings using the local client.
    """

    def getClient(self):
        return client.LocalClient(self.backend)


class PagingMixin(object):
    """
    Tests the paging code using a simulated backend.
    """
    @classmethod
    def setUpClass(cls):
        cls.numReferences = 25
        registry_db = utils.create_simulated_registry_db(
            random_seed=100, num_datasets=0,
            num_reference_sets=1,
            num_references_per_reference_set=cls.numReferences)
        cls.backend = backend.Backend(registry_db)
        cls.registry_db = registry_db

    def setUp(self):
        self.client = self.getClient()
        self.datamodelReferenceSet = self.registry_db.get_reference_sets()[0]
        self.datamodelReferences = self.datamodelReferenceSet.references
        self.references = [
            dmReference.get_protobuf()
            for dmReference in self.datamodelReferences]
        self.assertEqual(len(self.references), self.numReferences)

    def verifyAllReferences(self):
        """
        Verifies that we correctly return all references.
        """
        references = list(self.client.searchReferences(
            str(self.datamodelReferenceSet.id)))
        self.assertEqual(references, self.references)

    def testDefaultPageSize(self):
        self.verifyAllReferences()

    def verifyPageSize(self, pageSize):
        self.client.setPageSize(pageSize)
        self.assertEqual(pageSize, self.client.getPageSize())
        self.verifyAllReferences()

    def testPageSize1(self):
        self.verifyPageSize(1)

    def testPageSize2(self):
        self.verifyPageSize(2)

    def testPageSize3(self):
        self.verifyPageSize(3)

    def testPageSizeAlmostListLength(self):
        self.verifyPageSize(self.numReferences - 1)

    def testPageSizeListLength(self):
        self.verifyPageSize(self.numReferences)


class TestPagingLocal(PagingMixin, unittest.TestCase):
    """
    Tests paging using the local client.
    """

    def getClient(self):
        return client.LocalClient(self.backend)


class TestPagingHttp(PagingMixin, unittest.TestCase):
    """
    Tests paging using the HTTP client.
    """

    def getClient(self):
        return DummyHttpClient(self.backend)
