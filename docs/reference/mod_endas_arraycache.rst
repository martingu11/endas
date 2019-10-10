Large data cache (endas.arraycache)
===================================

Functionality for caching larger-than-memory data.

.. currentmodule:: endas.arraycache

Large array cache is used in situations where keeping all data in memory is likely not going to be feasible, such as in
Kalman Smoother implementations. Therefore, some data may be persisted to other storage (e.g. disk) for later retrieval.

The :class:`endas.arraycache.ArrayCache` class both defines the interface and serves as a trivial implementation.
:class:`ArrayCache` keeps all data in memory at all times and is therefore only suitable for when it is known that the
amount of data that needs to be held is going to be reasonably small. For more demanding situations, use other
implementation. It is also easy to write your own, if needed.

EnDAS comes with the following large array cache implementations:

.. autosummary::
    :toctree: generated/
    :template: class.rst
    :nosignatures:

    ArrayCache
    LRUArrayCache







