"""
Compiles the protobuf definitions.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import subprocess


def main():
    args = (
        "protoc -I . --python_out . "
        "ga4gh/proto/ga4gh/*.proto ga4gh/proto/google/protobuf/*.proto"
    )
    subprocess.check_call(args, shell=True)


if __name__ == "__main__":
    main()
