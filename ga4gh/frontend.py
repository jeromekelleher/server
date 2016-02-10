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


MIMETYPE = "application/json"
SEARCH_ENDPOINT_METHODS = ['POST', 'OPTIONS']
SECRET_KEY_LENGTH = 24

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
    return handler(body, acceptEncoding)

class Ga4ghProtocol(object):

    def __init__(self):
        repoDir = "ga4gh-example-data"
        repo = datarepo.FileSystemDataRepository(repoDir)
        self._backend = backend.Backend(repo)
        self.variantsets = VariantSet(self._backend)
        self.variants = Variant(self._backend)
        self.datasets = Dataset(self._backend)

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
