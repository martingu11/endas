#include <Endas/DA/IndexedPartitionPointQuery.hpp>
#include "../Compatibility.hpp"

#include <spatial/idle_point_multiset.hpp>
#include <spatial/region_iterator.hpp>

using namespace std;
using namespace endas;


struct compare_coord
{
    compare_coord(const Array2d& coords) : coordsRef(coords) { }
    bool operator() (spatial::dimension_type n, index_t a, index_t b) const
    {
        return coordsRef(n, a) < coordsRef(n, b);
    }
private:
    const Array2d& coordsRef;
};

struct coord_range_predicate
{
    coord_range_predicate(const AABox& query, const Array2d& coords) 
    : queryRef(query), coordsRef(coords) { }

    spatial::relative_order operator()(spatial::dimension_type n, spatial::dimension_type, 
                                       index_t x) const
    {
        double cn = coordsRef(n, x);
        double nmin = queryRef.min()(n);
        double nmax = queryRef.max()(n);
        return (cn < nmin)? spatial::below : ((cn >= nmax)? spatial::above : spatial::matching);
    }

private:
    const AABox& queryRef;
    const Array2d& coordsRef;  
};

typedef spatial::idle_point_multiset<0, index_t, compare_coord> PointIndex;



struct IndexedPartitionPointQuery::Data
{
    shared_ptr<const AABoxDomainPartitioning> partitioner;
    Array2d coords;

    PointIndex pointIndex; // Must be declared after `coords`!

    Data(shared_ptr<const AABoxDomainPartitioning> part, Array2d c)
    : partitioner(part), coords(move(c)),
      pointIndex(part->partitionCoordDim(), compare_coord(this->coords))
    { 
        ENDAS_ASSERT(partitioner);

        // The coordinate array must have the right number or spatial dimensions
        ENDAS_ASSERT(partitioner->partitionCoordDim() == coords.rows());

        // Insert all coords into the index. The index actually stores the indices and coordinates remain
        // in the given array. 
        // Todo: We may be getting lots of cache misses with this approach, if this is the case consider 
        // inserting the coordinates into the kd-tree

        auto n = coords.cols();
        vector<index_t> coordIndices;
        coordIndices.reserve(n);

        for (index_t i = 0; i != n; i++) coordIndices.push_back(i);
        pointIndex.insert_rebalance(coordIndices.begin(), coordIndices.end());
    }

};


IndexedPartitionPointQuery::IndexedPartitionPointQuery(shared_ptr<const AABoxDomainPartitioning> partitioner,
                                                       Array2d coords)
: mData(make_unique<Data>(partitioner, move(coords)))
{  }


IndexedPartitionPointQuery::~IndexedPartitionPointQuery()
{ }


void IndexedPartitionPointQuery::rangeQuery(int domain, double range, IndexArray& out) const
{
    // Find the bounding box of the domain
    AABox bbox = mData->partitioner->getLocalBox(domain);

    // Inflate the bounding box by the range
    bbox.min().array() -= range;
    bbox.max().array() += range;

    // Find all points inside the box. Use the exact distance from the box 
    double rangeSquared = range*range;
    coord_range_predicate predicate(bbox, mData->coords);

    for (auto iter = spatial::region_begin(mData->pointIndex, predicate);
         iter != spatial::region_end(mData->pointIndex, predicate); 
         ++iter)
    {
        index_t i = *iter;
        double distSquared = bbox.squaredExteriorDistance(mData->coords.col(i).matrix());

        if (distSquared <= rangeSquared) 
            out.push_back(i);
    }

}

