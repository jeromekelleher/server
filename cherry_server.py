"""
Experimental server based on cherrypy
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import cherrypy
from paste.translogger import TransLogger

import ga4gh.frontend as frontend


def run_server():
    cherrypy.config.update(
        {
            'server.socket_port': 8000,
            'server.socket_host': '0.0.0.0',
            'response.stream': True,
            # SSL stuff
            # 'server.ssl_module': 'builtin',
            'server.ssl_certificate': "cert.pem",
            'server.ssl_private_key': "privkey.pem",
    })
    cherrypy.quickstart(frontend.Ga4ghProtocol())

if __name__ == "__main__":
    run_server()



