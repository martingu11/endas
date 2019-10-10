"""
Functionality for caching larger-than-memory data.
"""

import os
import copy
import collections
import tempfile
import uuid

import numpy as np


class ArrayCache(object):
    """
    Trivial large data cache implementation.

    Large data cache is used in situations where keeping all data in memory is likely not going
    to be feasible. Therefore, some data may be persisted to other storage (e.g. disk) for later retrieval.

    Note:
      This implementation is basic and simply stores all data in memory via a Python dictionary. Also,
      it is not thread safe. Use one of subclasses for more suitable behaviour or implement your own.
    """

    def __init__(self):
        self._cache = {}
        self._keycounter = 0


    def put(self, data):
        """
        Places an object into the cache and returns handle for retrieval. The data is always copied.

        Args:
          data : Data instance to store.

        Returns : Handle object that represents the stored array.
        """
        self._keycounter+= 1
        self._cache[self._keycounter] = copy.deepcopy(data)
        return self._keycounter


    def get(self, handle, force_copy=False):
        """
        Retrieves data from the cache.

        Args:
          handle : Handle to the data instance to retrieve.

        Returns : Original data object.
        """
        x = self._cache[handle]
        if not force_copy: return x
        else: return copy.deepcopy(x)


    def get_exclusive(self, handle):
        """
        Retrieves data from the cache with exclusive and writable access.

        Args:
          handle : Handle to the data instance to retrieve.

        Returns : Original data object.
        """
        return self._cache[handle]


    def release(self, handle):
        """
        Releases exclusive access to object pointed by given handle.
        Only handles obtained by the `get_exclusive()` method should be released.

        Args:
          handle : Handle to the data instance to retrieve.

        Returns : None
        """
        pass


    def remove(self, handle):
        """
        Removes data from the cache.

        Args:
          arrayhandle : Handle to the data instance to remove.

        Returns: ``None``
        """
        del self._cache[handle]

    def clear(self):
        """
        Removes all data form the cache. This also invalidates all previously
        returned handles.
        """
        self._cache = {}
        self._keycounter = 0



class LRUArrayCache(ArrayCache):
    """
    Data cache that swaps recently used data objects to disk.

    The cache will keep items in memory as long as their combined size os below maxsizeMB. When the
    capacity is reached, least recently accessed items are swapped to storage to make space.

    This implementation uses file-based storage in the given `tempdir` or system default temp directory
    if `tempdir` is None. Subclasses can override retire(), restore() and drop() to implement other storage
    mechanisms.

    # Todo: Implement thread safety at some point.
    """

    def __init__(self, maxsizeMB=1024, tempdir=None):
        self._memoryitems = collections.OrderedDict()
        self._retireditems = {}

        self._keycounter = 0
        self._memorySizeMB = 0
        self._retiredSizeMB = 0
        self._capacityMB = maxsizeMB
        self._tempdir = tempdir if tempdir is not None else os.path.join(tempfile.gettempdir(), "dfscache")


        if not os.path.exists(self._tempdir):
            # For user-given temp directory require that it exists
            if tempdir is not None:
                raise RuntimeError("User-defined large array cache directory {} does not exist.".format(tempdir))
            # The default one is created automatically though
            else:
                os.makedirs(self._tempdir)

        #self._logger = config.getLogger("LRUDataCache")

    @property
    def tempdir(self): return self._tempdir


    def put(self, data):
        self._keycounter += 1
        datasize = self.itemsize(data)

        #self._logger.debug("Inserting data of size {}MB, assigned key {}".format(datasize, self._keycounter))


        # Need to make space for the new data item -> pop last-recently-used item(s) from the dictionary
        self._reservespace(datasize)

        self._memoryitems[self._keycounter] = (copy.deepcopy(data), datasize)
        self._memorySizeMB+= datasize

        return self._keycounter


    def get(self, handle, force_copy=False):
        # Item available in memory, just move it to the end of the memory list (most recently used)
        if handle in self._memoryitems:
            self._memoryitems.move_to_end(handle)
            data, size = next(reversed(self._memoryitems.values())) # Return last item
            return copy.deepcopy(data)

        # Not in memory? Perhaps it's been retired
        if handle in self._retireditems:
            retired_uuid, retired_size = self._retireditems[handle] # We'll keep it in storage also
            self._reservespace(retired_size)

            #self._logger.debug("Restoring item of size {}MB, key {}, uid {} from storage".format(retired_size, handle, retired_uuid))
            data = self.restore(retired_uuid)
            self._memoryitems[handle] = (data, retired_size)
            self._memorySizeMB+= retired_size
            return copy.deepcopy(data)

        raise KeyError("Item with handle {} not found in cache".format(handle))


    def get_exclusive(self, handle):
        raise NotImplementedError()

    def release(self, handle):
        raise NotImplementedError()


    def remove(self, handle):
        if handle in self._memoryitems:
            data, size = self._memoryitems.pop(handle)
            #self._logger.debug("Removing item of size {}MB, key {} from memory".format(size, handle))
            del data
            self._memorySizeMB-= size

        if handle in self._retireditems:
            retired_uuid, retired_size = self._retireditems.pop(handle)
            #self._logger.debug("Removing item of size {}MB, key {}, uuid {} from storage".format(retired_size, handle, retired_uuid))

            try: self.drop(retired_uuid)
            except: pass
            self._retiredSizeMB-= retired_size
            #self._logger.debug("Total occupied storage {}MB".format(self._retiredSizeMB))



    def clear(self):
        self._memoryitems.clear()
        self._memorySizeMB = 0

        for  retired_uuid, retired_size in self._retireditems.values():
            try: self.drop(retired_uuid)
            except: pass
        self._retireditems.clear()
        self._retiredSizeMB = 0


    def __enter__(self): return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear()

    def __del__(self):
        self.clear()


    def itemsize(self, data):
        if isinstance(data, np.ndarray): return data.nbytes // 1024 // 1024
        else: raise TypeError("LRUDataCache only supports numpy.ndarray instances.")

    def retire(self, data, uid):
        assert isinstance(data, np.ndarray)
        tmpfile = os.path.join(self._tempdir, str(uid)+'.npy')
        np.save(tmpfile, data, allow_pickle=False, fix_imports=False)


    def restore(self, uid):
        tmpfile = os.path.join(self._tempdir, str(uid)+'.npy')
        data = np.load(tmpfile, allow_pickle=False, fix_imports=False)
        assert isinstance(data, np.ndarray)
        return data


    def drop(self, uid):
        tmpfile = os.path.join(self._tempdir, str(uid) + '.npy')
        os.remove(tmpfile)


    def _reservespace(self, sizeMB):

        #self._logger.debug("Free memory space {}MB, need {}MB".format(
        #  max(self._capacityMB - self._memorySizeMB, 0),
        #  sizeMB))

        while self._memorySizeMB + sizeMB > self._capacityMB and len(self._memoryitems) > 0:
            lruitem_key, lruitem_value = self._memoryitems.popitem(last=False)
            lruitem_data, lruitem_size = lruitem_value

            if lruitem_key not in self._retireditems:
                retired_uuid = uuid.uuid4()
                #self._logger.debug("Retiring item of size {}MB, key {} to storage with uid {}".format(
                #  lruitem_size, lruitem_key, retired_uuid))

                self.retire(lruitem_data, retired_uuid)
                self._retireditems[lruitem_key] = (retired_uuid, lruitem_size)
                self._retiredSizeMB+= lruitem_size
                #self._logger.debug("Total occupied storage {}MB".format(self._retiredSizeMB))

            #else:
            #self._logger.debug("Retiring item of size {}MB, key {}, already in storage".format(
            #  lruitem_size, lruitem_key))

            del lruitem_data
            self._memorySizeMB -= lruitem_size


