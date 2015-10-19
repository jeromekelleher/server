"""
Command line interface programs for the GA4GH reference implementation.

TODO: document how to use these for development and simple deployment.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
import argparse
import sys
import logging
import unittest
import unittest.loader
import unittest.suite

import ga4gh.client as client
import ga4gh.protocol as protocol
import ga4gh.converters as converters
import ga4gh.frontend as frontend
import ga4gh.configtest as configtest


# the maximum value of a long type in avro = 2**63 - 1
# (64 bit signed integer)
# http://avro.apache.org/docs/1.7.7/spec.html#schema_primitive
# AVRO_LONG_MAX = (1 << 63) - 1
# TODO in the meantime, this is the max value pysam can handle
# This should be removed once pysam input sanitisation has been
# implemented.
AVRO_LONG_MAX = 2**31 - 1


def setCommaSeparatedAttribute(request, args, attr):
    # TODO document this function and possibly rename it as
    # part of the proto conversion.
    attribute = getattr(args, attr)
    if attribute is not None:
        proto_attr = getattr(request, attr)
        proto_attr.extend(attribute.split(","))


class RequestFactory(object):
    """
    Provides methods for easy initialization of request objects
    """
    def __init__(self, args):
        self.args = args
        self.workarounds = getWorkarounds(args)

    def usingWorkaroundsFor(self, workaround):
        """
        Returns true if we are using the passed-in workaround
        """
        return workaround in self.workarounds

    def createSearchVariantSetsRequest(self):
        request = protocol.SearchVariantSetsRequest()
        setCommaSeparatedAttribute(request, self.args, 'dataset_ids')
        request.page_size = self.args.page_size
        return request

    def createSearchVariantsRequest(self):
        request = protocol.SearchVariantsRequest()
        if self.args.reference_name is not None:
            request.reference_name = self.args.reference_name
        if self.args.variant_name is not None:
            request.variant_name = self.args.variant_name
        if self.args.start is not None:
            request.start = self.args.start
        if self.args.end is not None:
            request.end = self.args.end
        if self.usingWorkaroundsFor(client.HttpClient.workaroundGoogle):
            request.maxCalls = self.args.maxCalls
        if self.args.call_set_ids not in ['*', []]:
            request.call_set_ids = self.args.call_set_ids.split(",")
        setCommaSeparatedAttribute(request, self.args, 'variant_set_ids')
        return request

    def createSearchReferenceSetsRequest(self):
        request = protocol.SearchReferenceSetsRequest()
        setCommaSeparatedAttribute(request, self.args, 'accessions')
        setCommaSeparatedAttribute(request, self.args, 'md5checksums')
        return request

    def createSearchReferencesRequest(self):
        request = protocol.SearchReferencesRequest()
        setCommaSeparatedAttribute(request, self.args, 'accessions')
        setCommaSeparatedAttribute(request, self.args, 'md5checksums')
        return request

    def createSearchReadGroupSetsRequest(self):
        request = protocol.SearchReadGroupSetsRequest()
        setCommaSeparatedAttribute(request, self.args, 'dataset_ids')
        request.name = self.args.name
        return request

    def createSearchCallSetsRequest(self):
        request = protocol.SearchCallSetsRequest()
        setCommaSeparatedAttribute(request, self.args, 'variant_set_ids')
        request.name = self.args.name
        return request

    def createSearchReadsRequest(self):
        request = protocol.SearchReadsRequest()
        if self.usingWorkaroundsFor(client.HttpClient.workaroundGoogle):
            # google says reference_id not a valid field
            request = self.SearchReadsRequestGoogle()
        setCommaSeparatedAttribute(request, self.args, 'read_group_ids')
        request.start = self.args.start
        request.end = self.args.end
        request.reference_id = self.args.reference_id
        request.reference_name = self.args.reference_name
        return request

    def createSearchDatasetsRequest(self):
        request = protocol.SearchDatasetsRequest()
        return request

    def createListReferenceBasesRequest(self):
        request = protocol.ListReferenceBasesRequest()
        request.start = self.args.start
        request.end = self.args.end
        return request


def getWorkarounds(args):
    if args.workarounds is None:
        return set()
    else:
        return set(args.workarounds.split(','))


##############################################################################
# ga2vcf
##############################################################################


def addOutputFileArgument(parser):
    parser.add_argument(
        "--outputFile", "-o", default=None,
        help="the file to write the output to")


def ga2vcf_main(parser=None):
    # parse args
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH VCF conversion tool")
    addClientGlobalOptions(parser)
    addOutputFileArgument(parser)
    subparsers = parser.add_subparsers(title='subcommands',)
    variantsSearchParser = addVariantsSearchParser(subparsers)
    # TODO kinda ugly that we need to overload the parser with arguments
    # to populate the searchVariantSetsRequest;
    # this is a temporary workaround until we have getVariantSet
    # in the protocol;
    # also note both requests have a page_size attribute
    addDataSetIdsArgument(variantsSearchParser)
    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        ga2vcf_run(args)


def ga2vcf_run(args):
    # instantiate params
    searchVariantsRequest = RequestFactory(
        args).createSearchVariantsRequest()
    searchVariantSetsRequest = RequestFactory(
        args).createSearchVariantSetsRequest()
    if args.outputFile is None:
        outputStream = sys.stdout
    else:
        outputStream = open(args.outputFile, 'w')
    workarounds = getWorkarounds(args)
    httpClient = client.HttpClient(
        args.baseUrl, args.verbose, workarounds, args.key)

    # do conversion
    vcfConverter = converters.VcfConverter(
        httpClient, outputStream,
        searchVariantSetsRequest, searchVariantsRequest)
    vcfConverter.convert()

    # cleanup
    if args.outputFile is not None:
        outputStream.close()


##############################################################################
# ga2sam
##############################################################################


def ga2sam_main(parser=None):
    # parse args
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH SAM conversion tool")
    addClientGlobalOptions(parser)
    addReadsSearchParserArguments(parser)
    addBinaryOutputArgument(parser)
    addOutputFileArgument(parser)
    args = parser.parse_args()
    if "baseUrl" not in args:
        parser.print_help()
    else:
        ga2sam_run(args)


def ga2sam_run(args):
    # instantiate params
    searchReadsRequest = RequestFactory(
        args).createSearchReadsRequest()
    workarounds = getWorkarounds(args)
    httpClient = client.HttpClient(
        args.baseUrl, args.verbose, workarounds, args.key)

    # do conversion
    samConverter = converters.SamConverter(
        httpClient, searchReadsRequest, args.outputFile, args.binaryOutput)
    samConverter.convert()


def addBinaryOutputArgument(parser):
    parser.add_argument(
        "--binaryOutput", "-b", default=False,
        action="store_true", help="Output a BAM (binary) file")


##############################################################################
# Server
##############################################################################


def addGlobalOptions(parser):
    parser.add_argument(
        "--port", "-P", default=8000, type=int,
        help="The port to listen on")
    parser.add_argument(
        "--config", "-c", default='DevelopmentConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    parser.add_argument(
        "--dont-use-reloader", default=False, action="store_true",
        help="Don't use the flask reloader")


def server_main(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH reference server")
    addGlobalOptions(parser)
    args = parser.parse_args()
    frontend.configure(args.config_file, args.config)
    frontend.app.run(
        host="0.0.0.0", port=args.port,
        use_reloader=not args.dont_use_reloader)


##############################################################################
# Client
##############################################################################


class AbstractQueryRunner(object):
    """
    Abstract base class for runner classes
    """
    def __init__(self, args):
        self._workarounds = getWorkarounds(args)
        self._key = args.key
        self._verbosity = args.verbose
        self._httpClient = client.HttpClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)


class AbstractGetRunner(AbstractQueryRunner):
    """
    Abstract base class for get runner classes
    """
    def __init__(self, args):
        super(AbstractGetRunner, self).__init__(args)
        self._id = args.id
        self._httpClient = client.HttpClient(
            args.baseUrl, args.verbose, self._workarounds, self._key)

    def _run(self, method):
        response = method(self._id)
        print(response.id)


class AbstractSearchRunner(AbstractQueryRunner):
    """
    Abstract base class for search runner classes
    """
    def __init__(self, args):
        super(AbstractSearchRunner, self).__init__(args)

    def _setRequest(self, request, args):
        """
        Sets the _httpClient and other common attributes
        """
        self._minimalOutput = args.minimalOutput
        if 'page_size' in args:
            # ListReferenceBasesRequest does not have a page_size attr
            request.page_size = args.page_size
        self._request = request

    def _run(self, method, attrName=None):
        """
        Runs the request given methodname and prints out
        the each result's attrName attribute if it is provided.
        If not, prints each entire result object.
        """
        results = method(self._request)
        for result in results:
            if attrName is None:
                print(result)
            else:
                attr = getattr(result, attrName)
                print(attr)


class SearchVariantSetsRunner(AbstractSearchRunner):
    """
    Runner class for the variantsets/search method.
    """
    def __init__(self, args):
        super(SearchVariantSetsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchVariantSetsRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchVariantSets, 'id')


class SearchVariantsRunner(AbstractSearchRunner):
    """
    Runner class for the variants/search method.
    """
    def __init__(self, args):
        super(SearchVariantsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchVariantsRequest()
        # if no variantSets have been specified, send a request to
        # the server to grab all variantSets and then continue
        if args.variant_set_ids is None:
            datasetsRequest = protocol.SearchDatasetsRequest()
            datasetsResponse = self._httpClient.searchDatasets(
                datasetsRequest)
            dataset_ids = [dataset.id for dataset in datasetsResponse]
            variant_set_ids = []
            for dataset_id in dataset_ids:
                variantSetsRequest = protocol.SearchVariantSetsRequest()
                variantSetsRequest.dataset_ids = [dataset_id]
                response = self._httpClient.searchVariantSets(
                    variantSetsRequest)
                datasetVariantSetIds = [
                    variantSet.id for variantSet in response]
                variant_set_ids.extend(datasetVariantSetIds)
            request.variant_set_ids = variant_set_ids
        self._setRequest(request, args)

    def run(self):
        if self._minimalOutput:
            self._run(self._httpClient.searchVariants, 'id')
        else:
            results = self._httpClient.searchVariants(self._request)
            for result in results:
                self.printVariant(result)

    def printVariant(self, variant):
        """
        Prints out the specified Variant object in a VCF-like form.
        """
        print(
            variant.id, variant.variant_set_id, variant.names,
            variant.reference_name, variant.start, variant.end,
            variant.reference_bases, variant.alternate_bases,
            sep="\t", end="\t")
        for key, value in variant.info.items():
            print(key, value, sep="=", end=";")
        print("\t", end="")
        for c in variant.calls:
            info = "{"
            for k, v in c.info.items():
                info += "{}:{}, ".format(k, v)
            info += "}"
            print(
                c.genotype, c.genotype_likelihood, info,
                c.phaseset, sep=":", end="\t")
        print()


class SearchReferenceSetsRunner(AbstractSearchRunner):
    """
    Runner class for the referencesets/search method.
    """
    def __init__(self, args):
        super(SearchReferenceSetsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchReferenceSetsRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReferenceSets, 'id')


class SearchReferencesRunner(AbstractSearchRunner):
    """
    Runner class for the references/search method
    """
    def __init__(self, args):
        super(SearchReferencesRunner, self).__init__(args)
        request = RequestFactory(args).createSearchReferencesRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReferences, 'id')


class SearchReadGroupSetsRunner(AbstractSearchRunner):
    """
    Runner class for the readgroupsets/search method
    """
    def __init__(self, args):
        super(SearchReadGroupSetsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchReadGroupSetsRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReadGroupSets, 'id')


class SearchCallSetsRunner(AbstractSearchRunner):
    """
    Runner class for the callsets/search method
    """
    def __init__(self, args):
        super(SearchCallSetsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchCallSetsRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchCallSets, 'id')


class SearchReadsRunner(AbstractSearchRunner):
    """
    Runner class for the reads/search method
    """

    def __init__(self, args):
        super(SearchReadsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchReadsRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchReads, 'id')


class SearchDatasetsRunner(AbstractSearchRunner):
    """
    Runner class for the datasets/search method
    """
    def __init__(self, args):
        super(SearchDatasetsRunner, self).__init__(args)
        request = RequestFactory(args).createSearchDatasetsRequest()
        self._setRequest(request, args)

    def run(self):
        self._run(self._httpClient.searchDatasets, 'id')


class ListReferenceBasesRunner(AbstractSearchRunner):
    """
    Runner class for the references/{id}/bases method
    """
    def __init__(self, args):
        super(ListReferenceBasesRunner, self).__init__(args)
        request = RequestFactory(args).createListReferenceBasesRequest()
        self._id = args.id
        self._setRequest(request, args)

    def run(self):
        method = self._httpClient.listReferenceBases
        for base in method(self._request, self._id):
            print(base.sequence)


class GetReferenceSetRunner(AbstractGetRunner):
    """
    Runner class for the referencesets/{id} method
    """
    def __init__(self, args):
        super(GetReferenceSetRunner, self).__init__(args)

    def run(self):
        self._run(self._httpClient.getReferenceSet)


class GetReferenceRunner(AbstractGetRunner):
    """
    Runner class for the references/{id} method
    """
    def __init__(self, args):
        super(GetReferenceRunner, self).__init__(args)

    def run(self):
        self._run(self._httpClient.getReference)


class BenchmarkRunner(SearchVariantsRunner):
    """
    Runner class for the client side benchmarking. This is intended to give
    rough figures on protocol throughput on the server side over various
    requests.
    """
    def run(self):
        numVariants = 0
        beforeCpu = time.clock()
        beforeWall = time.time()
        try:
            for variant in self._httpClient.searchVariants(self._request):
                numVariants += 1
        except KeyboardInterrupt:
            pass
        cpuTime = time.clock() - beforeCpu
        wallTime = time.time() - beforeWall
        totalBytes = self._httpClient.getBytesRead()
        totalBytes /= 1024 * 1024
        s = "read {0} variants in {1:.2f} seconds; CPU time {2:.2f}".format(
            numVariants, wallTime, cpuTime)
        s += "; {0:.2f} MB @ {1:.2f} MB/s; {2:.2f} vars/s".format(
            totalBytes, totalBytes / wallTime, numVariants / wallTime)
        print(s)


def addVariantSearchOptions(parser):
    """
    Adds common options to a variant searches command line parser.
    """
    addVariantSetIdsArgument(parser)
    parser.add_argument(
        "--reference_name", "-r", default="1",
        help="Only return variants on this reference.")
    parser.add_argument(
        "--variant_name", "-n", default=None,
        help="Only return variants which have exactly this name.")
    parser.add_argument(
        "--call_set_ids", "-c", default=[],
        help="""Return variant calls which belong to call sets
            with these IDs. Pass in IDs as a comma separated list (no spaces).
            Omit this option to indicate 'all call sets'.
            """)
    addStartArgument(parser)
    addEndArgument(parser)
    addPageSizeArgument(parser)
    # maxCalls not in protocol; supported by google
    parser.add_argument(
        "--maxCalls", default=1,
        help="The maxiumum number of calls to return")


def addVariantSetIdsArgument(parser):
    parser.add_argument(
        "--variant_set_ids", "-V", default=None,
        help="The variant set id(s) to search over")


def addStartArgument(parser):
    parser.add_argument(
        "--start", "-s", default=0, type=int,
        help="The start of the search range (inclusive).")


def addEndArgument(parser):
    parser.add_argument(
        "--end", "-e", default=AVRO_LONG_MAX, type=int,
        help="The end of the search range (exclusive).")


def addIdArgument(parser):
    parser.add_argument("id", default=None, help="The id of the object")


def addGetArguments(parser):
    addIdArgument(parser)
    addUrlArgument(parser)


def addUrlArgument(parser):
    """
    Adds the URL endpoint argument to the specified parser.
    """
    parser.add_argument("baseUrl", help="The URL of the API endpoint")


def addAccessionsArgument(parser):
    parser.add_argument(
        "--accessions", default=None,
        help="The accessions to search over")


def addMd5ChecksumsArgument(parser):
    parser.add_argument(
        "--md5checksums", default=None,
        help="The md5checksums to search over")


def addPageSizeArgument(parser):
    parser.add_argument(
        "--page_size", "-m", default=100, type=int,
        help="The maximum number of results returned in one response.")


def addDataSetIdsArgument(parser):
    parser.add_argument(
        "--dataset_ids", default=None,
        help="The dataset_ids to search over")


def addNameArgument(parser):
    parser.add_argument(
        "--name", default=None,
        help="The name to search over")


def addClientGlobalOptions(parser):
    parser.add_argument('--verbose', '-v', action='count', default=0)
    parser.add_argument(
        "--workarounds", "-w", default=None, help="The workarounds to use")
    parser.add_argument(
        "--key", "-k", help="The auth key to use")
    parser.add_argument(
        "--minimalOutput", "-O", default=False,
        help="Use minimal output; default False",
        action='store_true')


def addHelpParser(subparsers):
    parser = subparsers.add_parser(
        "help", description="ga4gh_client help",
        help="show this help message and exit")
    return parser


def addBenchmarkingParser(subparsers):
    parser = subparsers.add_parser(
        "benchmark",
        description="Run simple benchmarks on the various methods",
        help="Benchmark server performance")
    parser.set_defaults(runner=BenchmarkRunner)
    addUrlArgument(parser)
    addVariantSearchOptions(parser)
    return parser


def addVariantsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "variants-search",
        description="Search for variants",
        help="Search for variants.")
    parser.set_defaults(runner=SearchVariantsRunner)
    addUrlArgument(parser)
    addVariantSearchOptions(parser)
    return parser


def addVariantSetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "variantsets-search",
        description="Search for variantSets",
        help="Search for variantSets.")
    parser.set_defaults(runner=SearchVariantSetsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addDataSetIdsArgument(parser)
    return parser


def addReferenceSetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "referencesets-search",
        description="Search for reference_sets",
        help="Search for reference_sets")
    parser.set_defaults(runner=SearchReferenceSetsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addAccessionsArgument(parser)
    addMd5ChecksumsArgument(parser)
    parser.add_argument(
        "--assemblyId",
        help="The assembly id to search over")
    return parser


def addReferencesSearchParser(subparsers):
    parser = subparsers.add_parser(
        "references-search",
        description="Search for references",
        help="Search for references")
    parser.set_defaults(runner=SearchReferencesRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addAccessionsArgument(parser)
    addMd5ChecksumsArgument(parser)
    return parser


def addReadGroupSetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "readgroupsets-search",
        description="Search for readGroupSets",
        help="Search for readGroupSets")
    parser.set_defaults(runner=SearchReadGroupSetsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addDataSetIdsArgument(parser)
    addNameArgument(parser)
    return parser


def addCallsetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "callsets-search",
        description="Search for call_sets",
        help="Search for call_sets")
    parser.set_defaults(runner=SearchCallSetsRunner)
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addNameArgument(parser)
    addVariantSetIdsArgument(parser)
    return parser


def addReadsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "reads-search",
        description="Search for reads",
        help="Search for reads")
    parser.set_defaults(runner=SearchReadsRunner)
    addReadsSearchParserArguments(parser)
    return parser


def addDatasetsSearchParser(subparsers):
    parser = subparsers.add_parser(
        "datasets-search",
        description="Search for datasets",
        help="Search for datasets")
    parser.set_defaults(runner=SearchDatasetsRunner)
    addUrlArgument(parser)
    return parser


def addReadsSearchParserArguments(parser):
    addUrlArgument(parser)
    addPageSizeArgument(parser)
    addStartArgument(parser)
    addEndArgument(parser)
    parser.add_argument(
        "--read_group_ids", default=None,
        help="The read_group_ids to search over")
    parser.add_argument(
        "--reference_id", default=None,
        help="The reference_id to search over")
    parser.add_argument(
        "--reference_name", default=None,
        help="The reference_name to search over")


def addReferenceSetsGetParser(subparsers):
    parser = subparsers.add_parser(
        "referencesets-get",
        description="Get a referenceset",
        help="Get a referenceset")
    parser.set_defaults(runner=GetReferenceSetRunner)
    addGetArguments(parser)


def addReferencesGetParser(subparsers):
    parser = subparsers.add_parser(
        "references-get",
        description="Get a reference",
        help="Get a reference")
    parser.set_defaults(runner=GetReferenceRunner)
    addGetArguments(parser)


def addReferencesBasesListParser(subparsers):
    parser = subparsers.add_parser(
        "references-list-bases",
        description="List bases of a reference",
        help="List bases of a reference")
    parser.set_defaults(runner=ListReferenceBasesRunner)
    addGetArguments(parser)
    addStartArgument(parser)
    addEndArgument(parser)


def client_main(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH reference client")
    addClientGlobalOptions(parser)
    subparsers = parser.add_subparsers(title='subcommands',)
    addHelpParser(subparsers)
    addBenchmarkingParser(subparsers)
    addVariantsSearchParser(subparsers)
    addVariantSetsSearchParser(subparsers)
    addReferenceSetsSearchParser(subparsers)
    addReferencesSearchParser(subparsers)
    addReadGroupSetsSearchParser(subparsers)
    addCallsetsSearchParser(subparsers)
    addReadsSearchParser(subparsers)
    addDatasetsSearchParser(subparsers)
    addReferenceSetsGetParser(subparsers)
    addReferencesGetParser(subparsers)
    addReferencesBasesListParser(subparsers)

    args = parser.parse_args()
    if "runner" not in args:
        parser.print_help()
    else:
        runner = args.runner(args)
        runner.run()

##############################################################################
# Configuration testing
##############################################################################


class SimplerResult(unittest.TestResult):
    """
    The TestResult class gives formatted tracebacks as error messages, which
    is not what we want. Instead we just want the error message from the
    err praram. Hence this subclass.
    """
    def addError(self, test, err):
        self.errors.append((test,
                            "{0}: {1}".format(err[0].__name__, err[1])))

    def addFailure(self, test, err):
        self.failures.append((test,
                              "{0}: {1}".format(err[0].__name__, err[1])))


def configtest_main(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="GA4GH server configuration validator")
    parser.add_argument(
        "--config", "-c", default='DevelopmentConfig', type=str,
        help="The configuration to use")
    parser.add_argument(
        "--config-file", "-f", type=str, default=None,
        help="The configuration file to use")
    args = parser.parse_args()
    configStr = 'ga4gh.serverconfig:{0}'.format(args.config)

    configtest.TestConfig.configStr = configStr
    configtest.TestConfig.configFile = args.config_file
    configtest.TestConfig.configEnv = "GA4GH_CONFIGURATION"

    loader = unittest.TestLoader()
    tests = loader.loadTestsFromModule(configtest)
    results = SimplerResult()
    tests.run(results)

    logging.basicConfig(level=logging.INFO)
    log = logging.getLogger(__name__)
    log.info('{0} Tests run. {1} errors, {2} failures, {3} skipped'.
             format(results.testsRun,
                    len(results.errors),
                    len(results.failures),
                    len(results.skipped)))
    for result in results.errors:
        if result is not None:
            log.critical('Error: {0}: {1}'.format(result[0].id(), result[1]))
    for result in results.failures:
        if result is not None:
            log.critical('Failure: {0}: {1}'.format(result[0].id(), result[1]))
    for result in results.skipped:
        if result is not None:
            log.info('Skipped: {0}: {1}'.format(result[0].id(), result[1]))
