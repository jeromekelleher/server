"""
Definitions of the GA4GH protocol types.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import google.protobuf.json_format as json_format

from ga4gh.proto.ga4gh.common_pb2 import *  # NOQA
from ga4gh.proto.ga4gh.variants_pb2 import *  # NOQA
from ga4gh.proto.ga4gh.reads_pb2 import *  # NOQA
from ga4gh.proto.ga4gh.datasets_pb2 import *  # NOQA
from ga4gh.proto.ga4gh.references_pb2 import *  # NOQA

# A map of response objects to the name of the attribute used to
# store the values returned.
_valueListNameMap = {
    SearchVariantSetsResponse: "variant_sets",
    SearchVariantsResponse: "variants",
}

version = "0.5.proto_experiment"


def getValueListName(protocolResponseClass):
    """
    Returns the name of the attribute in the specified protocol class
    that is used to hold the values in a search response.
    """
    return _valueListNameMap[protocolResponseClass]


def convertDatetime(t):
    """
    Converts the specified datetime object into its appropriate protocol
    value. This is the number of milliseconds from the epoch.
    """
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = t - epoch
    millis = delta.total_seconds() * 1000
    return int(millis)


class SearchResponseBuilder(object):
    """
    A class to allow sequential building of SearchResponse objects.
    """
    def __init__(self, responseClass, page_size, maxResponseLength):
        """
        Allocates a new SearchResponseBuilder for the specified
        responseClass, user-requested page_size and the system mandated
        maxResponseLength (in bytes). The maxResponseLength is an
        approximate limit on the overall length of the serialized
        response.
        """
        self._responseClass = responseClass
        self._page_size = page_size
        self._maxResponseLength = maxResponseLength
        self._numElements = 0
        self._next_page_token = None
        self._protoObject = responseClass()
        self._valueListName = getValueListName(responseClass)

    def getPageSize(self):
        """
        Returns the page size for this SearchResponseBuilder. This is the
        user-requested maximum size for the number of elements in the
        value list.
        """
        return self._page_size

    def getMaxResponseLength(self):
        """
        Returns the approximate maximum response length. More precisely,
        this is the total length (in bytes) of the concatenated JSON
        representations of the values in the value list after which
        we consider the buffer to be full.
        """
        return self._maxResponseLength

    def getNextPageToken(self):
        """
        Returns the value of the next_page_token for this
        SearchResponseBuilder.
        """
        return self._next_page_token

    def setNextPageToken(self, next_page_token):
        """
        Sets the next_page_token to the specified value.
        """
        self._next_page_token = next_page_token

    def addValue(self, protocolElement):
        """
        Appends the specified protocolElement to the value list for this
        response.
        """
        self._numElements += 1
        attr = getattr(self._protoObject, self._valueListName)
        obj = attr.add()
        obj.CopyFrom(protocolElement)

    def isFull(self):
        """
        Returns True if the response buffer is full, and False otherwise.
        The buffer is full if either (1) the number of items in the value
        list is >= pageSize or (2) the total length of the serialised
        elements in the page is >= maxResponseLength.
        """
        return (
            self._numElements >= self._page_size or
            self._protoObject.ByteSize() >= self._maxResponseLength)

    def getSerializedResponse(self):
        """
        Returns a string version of the SearchResponse that has
        been built by this SearchResponseBuilder.
        """
        if self._next_page_token is None:
            self._next_page_token = ""
        self._protoObject.next_page_token = self._next_page_token
        # s = self._protoObject.SerializeToString()
        s = json_format.MessageToJson(self._protoObject)
        return s
