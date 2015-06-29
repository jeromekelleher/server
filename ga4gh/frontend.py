"""
The Flask frontend for the GA4GH API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import datetime

import humanize
import cherrypy

import ga4gh
import ga4gh.backend as backend
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions

MIMETYPE = "application/json"
SEARCH_ENDPOINT_METHODS = ['POST', 'OPTIONS']


def http_methods_allowed(methods=['GET', 'HEAD']):
    method = cherrypy.request.method.upper()
    if method not in methods:
         cherrypy.response.headers['Allow'] = ", ".join(methods)
         raise cherrypy.HTTPError(405)

cherrypy.tools.allow = cherrypy.Tool('on_start_resource', http_methods_allowed)

def handlePostRequest(handler):
    cl = cherrypy.request.headers['Content-Length']
    body = cherrypy.request.body.read(int(cl))
    return handler(body)

class Ga4ghProtocol(object):

    def __init__(self):
        dataSource = "ga4gh-example-data"
        theBackend = backend.FileSystemBackend(dataSource)
        self._backend = theBackend
        self.variantsets = VariantSet(theBackend)
        self.variants = Variant(theBackend)
        self.datasets = Dataset(theBackend)

    @cherrypy.expose()
    def index(self):
        return "GA4GH API"

class Dataset(object):

    def __init__(self, backend):
        self._backend = backend

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.searchDatasets)


class VariantSet(object):

    def __init__(self, backend):
        self._backend = backend

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.searchVariantSets)

class Variant(object):

    def __init__(self, backend):
        self._backend = backend

    @cherrypy.expose()
    @cherrypy.tools.allow(methods=SEARCH_ENDPOINT_METHODS)
    def search(self):
        return handlePostRequest(self._backend.searchVariants)
