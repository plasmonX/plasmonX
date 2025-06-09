.. _algorithm:

algorithm section
=================

The ``algorithm`` section defines the numerical strategy used to solve the response equations.

The section is **optional**. plasmonX will decide the best algorithm based on the CPU/OMP parameters and the maximum memory available. 

Valid keywords
--------------

+--------------------------+----------------+--------------------------------------------------------+-----------------+
| Keyword                  | Type           | Allowed Values                                         |  Default        |
+==========================+================+========================================================+=================+
| ``method``               | string         | ``inversion``, ``iterative``, ``iterative on the fly`` |  --             |
+--------------------------+----------------+--------------------------------------------------------+-----------------+
| ``parallel execution``   | string         | ``matrix``, ``frequencies``                            |  --             |
+--------------------------+----------------+--------------------------------------------------------+-----------------+
| ``number of iterations`` | integer        | –                                                      | ``1``           |
+--------------------------+----------------+--------------------------------------------------------+-----------------+
| ``gmres dimension``      | integer        | –                                                      | ``1000``        |
+--------------------------+----------------+--------------------------------------------------------+-----------------+
| ``tolerance``            | float or string| –                                                      | ``1.0e-5``      |
+--------------------------+----------------+--------------------------------------------------------+-----------------+
| ``adaptive tuning``      | string         | ``yes``, ``no``                                        | ``yes``         |
+--------------------------+----------------+--------------------------------------------------------+-----------------+

Notes
-----

- ``method`` determines the solving strategy: 

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``inversion`` : inversion of the response matrix by LU factorization through LaPACK [default for small systems].
   - ``iterative`` : GMRES Krylov-based solver by storing the matrices on memory.
   - ``iterative on the fly`` : GMRES Krylov-based solver by on the fly matrix-vector products [default for large systems].
- ``parallel execution`` determines what to parallelize: 

   .. raw:: html

      <div style="margin: 0; padding: 0;"></div>

   - ``frequencies`` : the external cycle on frequencies is parallelized.
   - ``matrix`` : the construction of the matrices is parallelized.
- ``number of iterations`` determines the number of external iterations for `iterative` and `iterative on the fly` algorithm [not currently exploited]. 
- ``gmres dimension`` determines the number of microiterations for `iterative` and `iterative on the fly` algorithm [default: ``1000``]. 
- ``tolerance`` determines the convergence threshold (default: ``1.0e-5``). This value is scaled for ensure RMSE convergence.
- ``adaptive tuning`` determines whether to allow **plasmonX** to change the input parameters for best execution. 

Example
--------

.. code-block:: yaml

   algorithm:
     method: iterative
     parallel execution: matrix
     gmres dimension: 1000
     tolerance: 1e-6
     adaptive tuning: 'no'
