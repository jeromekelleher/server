"""
Tests for the repo manager tool
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import shutil
import tempfile
import unittest

import ga4gh.exceptions as exceptions
import ga4gh.registry as registry
import ga4gh.cli as cli
import tests.paths as paths


class TestGetNameFromPath(unittest.TestCase):
    """
    Tests the method for deriving the default name of objects from file
    paths.
    """
    def testError(self):
        self.assertRaises(ValueError, cli.getNameFromPath, "")

    def testLocalDirectory(self):
        self.assertEqual(cli.getNameFromPath("no_extension"), "no_extension")
        self.assertEqual(cli.getNameFromPath("x.y"), "x")
        self.assertEqual(cli.getNameFromPath("x.y.z"), "x")

    def testFullPaths(self):
        self.assertEqual(cli.getNameFromPath("/no_ext"), "no_ext")
        self.assertEqual(cli.getNameFromPath("/x.y"), "x")
        self.assertEqual(cli.getNameFromPath("/x.y.z"), "x")
        self.assertEqual(cli.getNameFromPath("/a/no_ext"), "no_ext")
        self.assertEqual(cli.getNameFromPath("/a/x.y"), "x")
        self.assertEqual(cli.getNameFromPath("/a/x.y.z"), "x")

    def testUrls(self):
        self.assertEqual(cli.getNameFromPath("file:///no_ext"), "no_ext")
        self.assertEqual(cli.getNameFromPath("http://example.com/x.y"), "x")
        self.assertEqual(cli.getNameFromPath("ftp://x.y.z"), "x")

    def testDirectoryName(self):
        self.assertEqual(cli.getNameFromPath("/a/xy"), "xy")
        self.assertEqual(cli.getNameFromPath("/a/xy/"), "xy")
        self.assertEqual(cli.getNameFromPath("xy/"), "xy")
        self.assertEqual(cli.getNameFromPath("xy"), "xy")


class AbstractRepoManagerTest(unittest.TestCase):
    """
    Base class for repo manager tests
    """
    def setUp(self):
        fd, path = tempfile.mkstemp(prefix="ga4gh_repoman_test")
        os.unlink(path)
        self._repoPath = "sqlite:///" + path
        self._repoFile = path

    def runCommand(self, cmd):
        cli.RepoManager.runCommand(cmd.split())

    def tearDown(self):
        os.unlink(self._repoFile)

    def readRepo(self):
        repo = registry.RegistryDb(self._repoPath)
        repo.open()
        return repo

    def init(self):
        self.runCommand("init {} -f".format(self._repoPath))

    def addOntology(self):
        # Add the sequence ontology
        self._ontologyName = paths.ontologyName
        cmd = "add-ontology {} {}".format(self._repoPath, paths.ontologyPath)
        self.runCommand(cmd)

    def addDataset(self, datasetName=None):
        if datasetName is None:
            datasetName = "test_dataset"
            self._datasetName = datasetName
        cmd = "add-dataset {} {}".format(self._repoPath, datasetName)
        self.runCommand(cmd)

    def addReferenceSet(self):
        self._referenceSetName = "test_rs"
        fastaFile = paths.faPath
        self.runCommand("add-referenceset {} {} --name={}".format(
            self._repoPath, fastaFile, self._referenceSetName))

    def addReadGroupSet(self):
        bamFile = paths.bamPath
        self._readGroupSetName = "test_rgs"
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, bamFile,
            self._referenceSetName, self._readGroupSetName)
        self.runCommand(cmd)

    def addVariantSet(self):
        vcfDir = paths.vcfDirPath
        self._variantSetName = "test_vs"
        cmd = (
            "add-variantset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, vcfDir,
            self._referenceSetName, self._variantSetName)
        self.runCommand(cmd)

    def addFeatureSet(self):
        featuresPath = paths.featuresPath
        self._featureSetName = paths.featureSetName
        cmd = (
            "add-featureset {} {} {} --referenceSetName={} "
            "--ontologyName={}").format(
            self._repoPath, self._datasetName, featuresPath,
            self._referenceSetName, self._ontologyName)
        self.runCommand(cmd)

    def getFeatureSet(self):
        repo = self.readRepo()
        dataset = repo.get_dataset_by_name(self._datasetName)
        featureSet = dataset.getFeatureSetByName(self._featureSetName)
        return featureSet


class TestAddFeatureSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddFeatureSet, self).setUp()
        self.init()
        self.addDataset()
        self.addOntology()
        self.addReferenceSet()

    def testAddFeatureSet(self):
        self.addFeatureSet()
        featureSet = self.getFeatureSet()
        self.assertEqual(featureSet.name, self._featureSetName)
        self.assertEqual(
            featureSet._parentContainer.name, self._datasetName)
        self.assertEqual(
            featureSet.getReferenceSet().name,
            self._referenceSetName)
        # TODO not clear these fields get populated now
        # self.assertEqual(featureSet.getInfo(), "TODO")
        # self.assertEqual(featureSet.getSourceUrl(), "TODO")

    def testAddFeatureSetNoReferenceSet(self):
        featuresPath = paths.featuresPath
        cmd = "add-featureset {} {} {} --ontologyName={}".format(
            self._repoPath, self._datasetName, featuresPath,
            self._ontologyName)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)

    def testAddFeatureSetBadReferenceSet(self):
        featuresPath = paths.featuresPath
        cmd = (
            "add-featureset {} {} {} --referenceSetName=notafefset "
            "--ontologyName={}").format(
            self._repoPath, self._datasetName, featuresPath,
            self._ontologyName)
        self.assertRaises(
            exceptions.ReferenceSetNameNotFoundException,
            self.runCommand, cmd)

    def testAddFeatureSetNoOntology(self):
        featuresPath = paths.featuresPath
        cmd = "add-featureset {} {} {} --referenceSetName={} ".format(
            self._repoPath, self._datasetName, featuresPath,
            self._referenceSetName)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)

    def testAddFeatureSetBadOntology(self):
        featuresPath = paths.featuresPath
        cmd = "add-featureset {} {} {} --referenceSetName={} ".format(
            self._repoPath, self._datasetName, featuresPath,
            self._referenceSetName)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)


class TestRemoveFeatureSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveFeatureSet, self).setUp()
        self.init()
        self.addDataset()
        self.addOntology()
        self.addReferenceSet()
        self.addFeatureSet()

    def testRemoveFeatureSet(self):
        featureSet = self.getFeatureSet()
        cmd = "remove-featureset {} {} {} -f".format(
            self._repoPath, self._datasetName, featureSet.name)
        self.runCommand(cmd)
        with self.assertRaises(exceptions.FeatureSetNameNotFoundException):
            self.getFeatureSet()


class TestAddDataset(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddDataset, self).setUp()
        self.init()

    def testDefaults(self):
        name = "test_dataset"
        self.runCommand("add-dataset {} {}".format(self._repoPath, name))
        repo = self.readRepo()
        dataset = repo.get_dataset_by_name(name)
        self.assertEqual(dataset.name, name)

    def testSameName(self):
        name = "test_dataset"
        cmd = "add-dataset {} {}".format(self._repoPath, name)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)


class TestAddReferenceSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddReferenceSet, self).setUp()
        self.init()

    def testDefaults(self):
        fasta_file = paths.ncbi37FaPath
        name = os.path.split(fasta_file)[1].split(".")[0]
        self.runCommand("add-referenceset {} {}".format(
            self._repoPath, fasta_file))
        repo = self.readRepo()
        reference_set = repo.get_reference_set_by_name(name)
        self.assertEqual(reference_set.name, name)
        for reference in reference_set.references:
            self.assertEqual(reference.fasta_url, os.path.abspath(fasta_file))
        # TODO check that the default values for all fields are set correctly.

    def testWithName(self):
        name = "test_reference_set"
        fasta_file = paths.ncbi37FaPath
        cmd = "add-referenceset {} {} --name={}".format(
            self._repoPath, fasta_file, name)
        self.runCommand(cmd)
        repo = self.readRepo()
        reference_set = repo.get_reference_set_by_name(name)
        self.assertEqual(reference_set.name, name)
        for reference in reference_set.references:
            self.assertEqual(reference.fasta_url, os.path.abspath(fasta_file))

    def testWithSameName(self):
        fastaFile = paths.ncbi37FaPath
        # Default name
        cmd = "add-referenceset {} {}".format(self._repoPath, fastaFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)
        # Specified name
        cmd = "add-referenceset {} {} --name=testname".format(
            self._repoPath, fastaFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)


class TestAddOntology(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddOntology, self).setUp()
        self.init()

    def testDefaults(self):
        ontologyFile = paths.ontologyPath
        name = os.path.split(ontologyFile)[1].split(".")[0]
        self.runCommand("add-ontology {} {}".format(
            self._repoPath, ontologyFile))
        repo = self.readRepo()
        ontology = repo.getOntologyByName(name)
        self.assertEqual(ontology.getName(), name)
        self.assertEqual(ontology.getDataUrl(), os.path.abspath(ontologyFile))

    def testWithName(self):
        ontologyFile = paths.ontologyPath
        name = "test_name"
        self.runCommand("add-ontology {} {} --name={}".format(
            self._repoPath, ontologyFile, name))
        repo = self.readRepo()
        ontology = repo.getOntologyByName(name)
        self.assertEqual(ontology.getName(), name)
        self.assertEqual(ontology.getDataUrl(), os.path.abspath(ontologyFile))

    def testWithSameName(self):
        ontologyFile = paths.ontologyPath
        # Default name
        cmd = "add-ontology {} {}".format(self._repoPath, ontologyFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)
        # Specified name
        cmd = "add-ontology {} {} --name=testname".format(
            self._repoPath, ontologyFile)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)

    def testMissingFile(self):
        cmd = "add-ontology {} {}".format(self._repoPath, "/no/such/file")
        self.assertRaises(
            exceptions.FileOpenFailedException, self.runCommand, cmd)

    def testNonOboTextFile(self):
        cmd = "add-ontology {} {}".format(
            self._repoPath, paths.landingMessageHtml)
        self.assertRaises(
            exceptions.OntologyFileFormatException, self.runCommand, cmd)

    def testNonOboBinaryFile(self):
        cmd = "add-ontology {} {}".format(self._repoPath, paths.bamPath)
        self.assertRaises(
            exceptions.OntologyFileFormatException, self.runCommand, cmd)


class TestRemoveDataset(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveDataset, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()

    def assertDatasetRemoved(self):
        repo = self.readRepo()
        self.assertRaises(
            exceptions.DatasetNameNotFoundException,
            repo.get_dataset_by_name, self._datasetName)

    def testEmptyDatasetForce(self):
        self.runCommand("remove-dataset {} {} -f".format(
            self._repoPath, self._datasetName))
        self.assertDatasetRemoved()

    def testContainsReadGroupSet(self):
        self.addReadGroupSet()
        self.runCommand("remove-dataset {} {} -f".format(
            self._repoPath, self._datasetName))
        self.assertDatasetRemoved()


class TestRemoveReadGroupSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveReadGroupSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()
        self.addReadGroupSet()

    def assertReadGroupSetRemoved(self):
        repo = self.readRepo()
        dataset = repo.get_dataset_by_name(self._datasetName)
        self.assertRaises(
            exceptions.ReadGroupSetNameNotFoundException,
            dataset.getReadGroupSetByName, self._readGroupSetName)

    def testWithForce(self):
        self.runCommand("remove-readgroupset {} {} {} -f".format(
            self._repoPath, self._datasetName, self._readGroupSetName))
        self.assertReadGroupSetRemoved()


class TestRemoveVariantSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveVariantSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()
        self.addVariantSet()

    def assertVariantSetRemoved(self):
        repo = self.readRepo()
        dataset = repo.get_dataset_by_name(self._datasetName)
        self.assertRaises(
            exceptions.VariantSetNameNotFoundException,
            dataset.get_variant_set_by_name, self._variantSetName)

    def testWithForce(self):
        self.runCommand("remove-variantset {} {} {} -f".format(
            self._repoPath, self._datasetName, self._variantSetName))
        self.assertVariantSetRemoved()

    # TODO test when we have a variant set with the same name in
    # a different dataset. This should be unaffected.


class TestRemoveReferenceSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveReferenceSet, self).setUp()
        self.init()
        self.addReferenceSet()

    def assertReferenceSetRemoved(self):
        repo = self.readRepo()
        self.assertRaises(
            exceptions.ReferenceSetNameNotFoundException,
            repo.get_reference_set_by_name, self._referenceSetName)

    def testDefaults(self):
        self.runCommand("remove-referenceset {} {} -f".format(
            self._repoPath, self._referenceSetName))
        self.assertReferenceSetRemoved()


class TestVerify(AbstractRepoManagerTest):

    def setUp(self):
        super(TestVerify, self).setUp()

    def testVerify(self):
        self.init()
        self.addDataset()
        self.addOntology()
        self.addReferenceSet()
        self.addReadGroupSet()
        self.addFeatureSet()
        self.addVariantSet()
        cmd = "verify {}".format(self._repoPath)
        self.runCommand(cmd)


class TestRemoveOntology(AbstractRepoManagerTest):

    def setUp(self):
        super(TestRemoveOntology, self).setUp()
        self.init()
        self.addOntology()

    def assertOntologyRemoved(self):
        repo = self.readRepo()
        self.assertRaises(
            exceptions.OntologyNameNotFoundException,
            repo.getOntologyByName, self._ontologyName)

    def testDefaults(self):
        self.runCommand("remove-ontology {} {} -f".format(
            self._repoPath, self._ontologyName))
        self.assertOntologyRemoved()


class TestAddReadGroupSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddReadGroupSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()

    def verifyReadGroupSet(self, name, dataUrl, indexFile):
        repo = self.readRepo()
        referenceSet = repo.get_reference_set_by_name(self._referenceSetName)
        readGroupSet = repo.get_read_group_set_by_name(self._datasetName, name)
        self.assertEqual(readGroupSet.name, name)
        self.assertEqual(readGroupSet.reference_set, referenceSet)
        self.assertEqual(readGroupSet.data_url, os.path.abspath(dataUrl))
        self.assertEqual(
            readGroupSet.index_url, os.path.abspath(indexFile))

    def testDefaultsLocalFile(self):
        bamFile = paths.bamPath
        name = os.path.split(bamFile)[1].split(".")[0]
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                self._referenceSetName)
        self.runCommand(cmd)
        self.verifyReadGroupSet(name, bamFile, bamFile + ".bai")

    def testLocalFileWithIndex(self):
        bamFile = paths.bamPath
        name = os.path.split(bamFile)[1].split(".")[0]
        with tempfile.NamedTemporaryFile() as temp:
            indexFile = temp.name
            shutil.copyfile(bamFile + ".bai", indexFile)
            cmd = (
                "add-readgroupset {} {} {} -I {} "
                "--referenceSetName={}").format(
                    self._repoPath, self._datasetName, bamFile,
                    indexFile, self._referenceSetName)
            self.runCommand(cmd)
            self.verifyReadGroupSet(name, bamFile, indexFile)

    def testLocalFileWithName(self):
        bamFile = paths.bamPath
        name = "test_rgs"
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, bamFile,
            self._referenceSetName, name)
        self.runCommand(cmd)
        self.verifyReadGroupSet(name, bamFile, bamFile + ".bai")

    def testAddReadGroupSetWithSameName(self):
        # Default name
        bamFile = paths.bamPath
        name = os.path.split(bamFile)[1].split(".")[0]
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                self._referenceSetName)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)
        # Specified name
        name = "test_rgs"
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, bamFile,
            self._referenceSetName, name)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)

    def testUrlWithMissingIndex(self):
        bamFile = "http://example.com/example.bam"
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                self._referenceSetName)
        self.assertRaises(
            exceptions.MissingIndexException, self.runCommand, cmd)

    def testMissingDataset(self):
        bamFile = paths.bamPath
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, "not_a_dataset_name", bamFile,
                self._referenceSetName)
        self.assertRaises(
            exceptions.DatasetNameNotFoundException, self.runCommand, cmd)

    def testMissingReferenceSet(self):
        bamFile = paths.bamPath
        cmd = "add-readgroupset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, bamFile,
                "not_a_referenceset_name")
        self.assertRaises(
            exceptions.ReferenceSetNameNotFoundException, self.runCommand, cmd)


class TestAddVariantSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddVariantSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()
        self.vcfDir = paths.vcfDirPath
        self.vcfFiles = glob.glob(os.path.join(paths.vcfDirPath, "*.vcf.gz"))
        self.indexFiles = [vcfFile + ".tbi" for vcfFile in self.vcfFiles]

    def verifyVariantSet(self, name, data_urls, index_files):
        repo = self.readRepo()
        reference_set = repo.get_reference_set_by_name(self._referenceSetName)
        variant_set = repo.get_variant_set_by_name(self._datasetName, name)
        self.assertEqual(variant_set.name, name)
        self.assertEqual(variant_set.reference_set, reference_set)
        data_urls = map(lambda x: os.path.abspath(x), data_urls)
        index_files = map(lambda x: os.path.abspath(x), index_files)
        pairs = sorted(zip(data_urls, index_files))
        other_pairs = sorted(
            (vcf_file.data_url, vcf_file.index_url)
            for vcf_file in variant_set.vcf_files)
        self.assertEqual(pairs, other_pairs)

    def testDefaultsLocalFiles(self):
        dataFiles = self.vcfFiles
        name = "test_name"
        cmd = "add-variantset {} {} {} --name={} --referenceSetName={}".format(
                self._repoPath, self._datasetName, " ".join(dataFiles),
                name, self._referenceSetName)
        self.runCommand(cmd)
        self.verifyVariantSet(name, dataFiles, self.indexFiles)

    def testDefaultsLocalDirectory(self):
        vcfDir = self.vcfDir
        name = os.path.split(vcfDir)[1]
        cmd = "add-variantset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, vcfDir,
                self._referenceSetName)
        self.runCommand(cmd)
        self.verifyVariantSet(name, self.vcfFiles, self.indexFiles)

    def testLocalFilesWithIndexes(self):
        dataFiles = self.vcfFiles
        tempdir = tempfile.mkdtemp(prefix="ga4gh_test_add_variantset")
        name = "test_name"
        try:
            indexFiles = []
            for indexFile in self.indexFiles:
                indexFileCopy = os.path.join(
                    tempdir, os.path.split(indexFile)[1])
                shutil.copyfile(indexFile, indexFileCopy)
                indexFiles.append(indexFileCopy)
            cmd = (
                "add-variantset {} {} {} -I {} --name={} "
                "--referenceSetName={}".format(
                    self._repoPath, self._datasetName, " ".join(dataFiles),
                    " ".join(indexFiles), name, self._referenceSetName))
            self.runCommand(cmd)
            self.verifyVariantSet(name, dataFiles, indexFiles)
        finally:
            shutil.rmtree(tempdir)

    def testAddVariantSetWithSameName(self):
        # Default name
        vcfDir = self.vcfDir
        name = os.path.split(vcfDir)[1]
        cmd = "add-variantset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, vcfDir,
                self._referenceSetName)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)
        # Specified name
        name = "test_vs"
        cmd = (
            "add-variantset {} {} {} --referenceSetName={} "
            "--name={}").format(
            self._repoPath, self._datasetName, vcfDir,
            self._referenceSetName, name)
        self.runCommand(cmd)
        self.assertRaises(
            exceptions.DuplicateNameException, self.runCommand, cmd)

    def testUrlWithMissingIndex(self):
        dataFile = "http://example.com/example.vcf.gz"
        cmd = "add-variantset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, dataFile,
                self._referenceSetName)
        self.assertRaises(
            exceptions.MissingIndexException, self.runCommand, cmd)

    def testMissingDataset(self):
        cmd = "add-variantset {} {} {} --referenceSetName={}".format(
                self._repoPath, "not_a_dataset_name", self.vcfDir,
                self._referenceSetName)
        self.assertRaises(
            exceptions.DatasetNameNotFoundException, self.runCommand, cmd)

    def testMissingReferenceSet(self):
        cmd = "add-variantset {} {} {} --referenceSetName={}".format(
                self._repoPath, self._datasetName, self.vcfDir,
                "not_a_referenceset_name")
        self.assertRaises(
            exceptions.ReferenceSetNameNotFoundException, self.runCommand, cmd)

    # TODO add more tests for to verify that errors are correctly thrown
    # when incorrect indexes are passed, mixtures of directories and URLS
    # for the dataFiles argument, and other common error cases in the UI.


class TestAddAnnotatedVariantSet(AbstractRepoManagerTest):

    def setUp(self):
        super(TestAddAnnotatedVariantSet, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()
        self.addOntology()
        self.vcfDir = paths.annotatedVcfPath

    def testNoAnnotations(self):
        name = "test_vs_no_annotations"
        cmd = "add-variantset {} {} {} -R {} -n {}".format(
            self._repoPath, self._datasetName, self.vcfDir,
            self._referenceSetName, name)
        self.runCommand(cmd)
        repo = self.readRepo()
        dataset = repo.get_dataset_by_name(self._datasetName)
        variantSet = dataset.get_variant_set_by_name(name)
        self.assertEqual(len(variantSet.getVariantAnnotationSets()), 0)

    def testAnnotations(self):
        name = "test_vs_annotations"
        cmd = "add-variantset {} {} {} -R {} -n {} -aO {}".format(
            self._repoPath, self._datasetName, self.vcfDir,
            self._referenceSetName, name, self._ontologyName)
        self.runCommand(cmd)
        repo = self.readRepo()
        dataset = repo.get_dataset_by_name(self._datasetName)
        variantSet = dataset.get_variant_set_by_name(name)
        self.assertEqual(len(variantSet.getVariantAnnotationSets()), 1)

    def testAnnotationsNoOntology(self):
        name = "test_vs_annotations"
        cmd = "add-variantset {} {} {} -R {} -n {} -a".format(
            self._repoPath, self._datasetName, self.vcfDir,
            self._referenceSetName, name)
        self.assertRaises(
            exceptions.RepoManagerException, self.runCommand, cmd)

    def testAnnotationsBadOntology(self):
        name = "test_vs_annotations"
        cmd = "add-variantset {} {} {} -R {} -n {} -aO {}".format(
            self._repoPath, self._datasetName, self.vcfDir,
            self._referenceSetName, name, "not_an_ontology")
        self.assertRaises(
            exceptions.OntologyNameNotFoundException, self.runCommand, cmd)


class TestDuplicateNameDelete(AbstractRepoManagerTest):
    """
    If two objects exist with the same name in different datasets,
    ensure that only one is deleted on a delete call
    """
    def setUp(self):
        super(TestDuplicateNameDelete, self).setUp()
        self.init()
        self.dataset1Name = "dataset1"
        self.dataset2Name = "dataset2"
        self.addDataset(self.dataset1Name)
        self.addDataset(self.dataset2Name)
        self.addOntology()
        self.addReferenceSet()

    def readDatasets(self):
        repo = self.readRepo()
        self.dataset1 = repo.get_dataset_by_name(self.dataset1Name)
        self.dataset2 = repo.get_dataset_by_name(self.dataset2Name)

    def testReadGroupSetDelete(self):
        readGroupSetName = "test_rgs"
        cmdString = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "--name={}")
        addReadGroupSetCmd1 = cmdString.format(
            self._repoPath, self.dataset1Name, paths.bamPath,
            self._referenceSetName, readGroupSetName)
        self.runCommand(addReadGroupSetCmd1)
        addReadGroupSetCmd2 = cmdString.format(
            self._repoPath, self.dataset2Name, paths.bamPath,
            self._referenceSetName, readGroupSetName)
        self.runCommand(addReadGroupSetCmd2)
        removeCmd = "remove-readgroupset {} {} {} -f".format(
            self._repoPath, self.dataset1Name, readGroupSetName)
        self.runCommand(removeCmd)
        self.readDatasets()
        self.assertEqual(len(self.dataset1.getReadGroupSets()), 0)
        self.assertEqual(len(self.dataset2.getReadGroupSets()), 1)

    def testVariantSetDelete(self):
        vcfDir = paths.vcfDirPath
        variantSetName = "test_vs"
        cmdString = "add-variantset {} {} {} --referenceSetName={} --name={}"
        addVariantSetCmd1 = cmdString.format(
            self._repoPath, self.dataset1Name, vcfDir,
            self._referenceSetName, variantSetName)
        self.runCommand(addVariantSetCmd1)
        addVariantSetCmd2 = cmdString.format(
            self._repoPath, self.dataset2Name, vcfDir,
            self._referenceSetName, variantSetName)
        self.runCommand(addVariantSetCmd2)
        removeCmd = "remove-variantset {} {} {} -f".format(
            self._repoPath, self.dataset1Name, variantSetName)
        self.runCommand(removeCmd)
        self.readDatasets()
        self.assertEqual(len(self.dataset1.getVariantSets()), 0)
        self.assertEqual(len(self.dataset2.getVariantSets()), 1)

    def testFeatureSetDelete(self):
        cmdString = "add-featureset {} {} {} -R {} -O {}"
        addFeatureSetCmd1 = cmdString.format(
            self._repoPath, self.dataset1Name, paths.featuresPath,
            self._referenceSetName, self._ontologyName)
        self.runCommand(addFeatureSetCmd1)
        addFeatureSetCmd2 = cmdString.format(
            self._repoPath, self.dataset2Name, paths.featuresPath,
            self._referenceSetName, self._ontologyName)
        self.runCommand(addFeatureSetCmd2)
        removeCmd = "remove-featureset {} {} {} -f".format(
            self._repoPath, self.dataset1Name, paths.featureSetName)
        self.runCommand(removeCmd)
        self.readDatasets()
        self.assertEqual(len(self.dataset1.getFeatureSets()), 0)
        self.assertEqual(len(self.dataset2.getFeatureSets()), 1)


class TestInvalidVariantIndexFile(AbstractRepoManagerTest):
    """
    Test that the repo manager throws exceptions when invalid index
    files are provided for vcf files.
    """
    def setUp(self):
        super(TestInvalidVariantIndexFile, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()

    def _testWithIndexPath(self, indexPath):
        cmd = (
            "add-variantset {} {} {} --referenceSetName={} -I {}").format(
                self._repoPath, self._datasetName, paths.vcfPath1,
                self._referenceSetName, indexPath)
        with self.assertRaises(exceptions.NotIndexedException):
            self.runCommand(cmd)

    def testNonexistentIndexFile(self):
        indexPath = '/path/does/not/exist'
        self._testWithIndexPath(indexPath)

    def testIndexFileNotAnIndexFile(self):
        indexPath = paths.vcfPath2  # not an index file
        self._testWithIndexPath(indexPath)

    @unittest.skip("Skipping until we can detect incorrect indexes")
    def testWrongIndexFile(self):
        indexPath = paths.vcfIndexPath2  # incorrect index
        self._testWithIndexPath(indexPath)


class TestInvalidReadGroupSetIndexFile(AbstractRepoManagerTest):
    """
    Test that the repo manager throws exceptions when invalid index
    files are provided for bam files.
    """
    @classmethod
    def setUpClass(cls):
        # clear the file handle cache because if the data file we are
        # testing with an invalid index is already in the cache, the
        # index will not be opened during the test; without this line
        # the below tests will succeed when the test class is run but
        # fail when the file's tests are run
        datamodel.fileHandleCache = datamodel.PysamFileHandleCache()

    def setUp(self):
        super(TestInvalidReadGroupSetIndexFile, self).setUp()
        self.init()
        self.addDataset()
        self.addReferenceSet()

    def _testWithIndexPath(self, indexPath):
        cmd = (
            "add-readgroupset {} {} {} --referenceSetName={} "
            "-I {}").format(
                self._repoPath, self._datasetName, paths.bamPath,
                self._referenceSetName, indexPath)
        self.runCommand(cmd)

    def testNonexistentIndexFile(self):
        indexPath = '/path/does/not/exist'
        with self.assertRaises(exceptions.FileOpenFailedException):
            self._testWithIndexPath(indexPath)

    def testIndexFileNotAnIndexFile(self):
        indexPath = paths.bamPath2  # not an index file
        with self.assertRaises(exceptions.DataException):
            self._testWithIndexPath(indexPath)

    @unittest.skip("Skipping until we can detect incorrect indexes")
    def testWrongIndexFile(self):
        indexPath = paths.bamIndexPath2  # incorrect index
        self._testWithIndexPath(indexPath)
