# -*- coding: utf-8 -*-
# Copyright (C) 2012, Robert Schroll
#
# Visvis is distributed under the terms of the (new) BSD License.
# The full license can be found in 'license.txt'.

import visvis as vv


if __name__ == '__main__':
    # Create BaseMesh object (has no visualization props)
    bm = vv.meshRead('teapot.ssdf')
    # Show it, returning a Mesh object (which does have visualization props)
    m = vv.mesh(bm)
    app = vv.use()
    app.Run()
