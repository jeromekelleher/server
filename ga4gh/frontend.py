"""
The Flask frontend for the GA4GH API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import datetime
import socket
import urlparse
import functools

import humanize
import cherrypy

import ga4gh
import ga4gh.backend as backend
import ga4gh.datamodel as datamodel
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions
import ga4gh.datarepo as datarepo


SEARCH_ENDPOINT_METHODS = ['POST', 'OPTIONS']
GET_ENDPOINT_METHODS = ['GET']

def http_methods_allowed(methods=['GET', 'HEAD']):
    method = cherrypy.request.method.upper()
    if method not in methods:
         cherrypy.response.headers['Allow'] = ", ".join(methods)
         raise cherrypy.HTTPError(405)

cherrypy.tools.allow = cherrypy.Tool('on_start_resource', http_methods_allowed)


def handlePostRequest(handler):
    cl = cherrypy.request.headers['Content-Length']
    body = cherrypy.request.body.read(int(cl))
    acceptEncoding = cherrypy.request.headers['Accept-Encoding']
    cherrypy.response.headers['Content-Type']= acceptEncoding
    return handler(body, acceptEncoding)

def handleGetRequest(handler, id_, key):
    # TODO handle niceties of HTTP stuff here.
    return handler(id_)


class Ga4ghProtocol(object):

    def __init__(self):
        repoDir = "ga4gh-example-data"
        repo = datarepo.FileSystemDataRepository(repoDir)
        self._backend = backend.Backend(repo)
        self.datasets = Dataset(self._backend)
        self.referencesets = ReferenceSet(self._backend)
        self.references = Reference(self._backend)
        self.variantsets = VariantSet(self._backend)
        self.variants = Variant(self._backend)
        self.readgroupsets = ReadGroupSet(self._backend)
        self.readgroups = ReadGroup(self._backend)
        self.reads = Read(self._backend)

    @cherrypy.expose()
    def index(self):
        return "GA4GH API"

class ServerObject(object):

    def __init__(self, backend):
        self._backend = backend

class Dataset(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchDatasets)

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=GET_ENDPOINT_METHODS)
    def default(self, id_, key=None):
        return handleGetRequest(self._backend.runGetDataset, id_, key)


# TODO add in support for other objects by
# 1) adding in a new class here;
# 2) updating backend function to use new pattern.

class ReferenceSet(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchReferenceSets)

class Reference(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchReferences)


class VariantSet(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchVariantSets)

class Variant(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchVariants)

class ReadGroupSet(ServerObject):
    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchReadGroupSets)

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=GET_ENDPOINT_METHODS)
    def default(self, id_, key=None):
        return handleGetRequest(self._backend.runGetReadGroupSet, id_, key)

class ReadGroup(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=GET_ENDPOINT_METHODS)
    def default(self, id_, key=None):
        return handleGetRequest(self._backend.runGetReadGroup, id_, key)



class Read(ServerObject):

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.runSearchReads)


