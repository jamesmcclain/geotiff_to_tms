#ifndef __RTREE_H__
#define __RTREE_H__

#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bi  = boost::interprocess;
namespace bg  = boost::geometry;
namespace bgm = boost::geometry::model;
namespace bgi = boost::geometry::index;

// Resource: http://www.boost.org/doc/libs/1_63_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/index_stored_in_mapped_file_using_boost_interprocess.html
typedef bgm::point<double, 2, bg::cs::cartesian> point_t;
typedef bgm::box<point_t> box_t;
typedef std::pair<box_t, lesser_landsat_scene_struct> value_t;
typedef bgi::linear<32,8> params_t;
typedef bgi::indexable<value_t> indexable_t;
typedef bgi::equal_to<value_t> equal_to_t;
typedef bi::allocator<value_t, bi::managed_mapped_file::segment_manager> allocator_t;
typedef bgi::rtree<value_t, params_t, indexable_t, equal_to_t, allocator_t> rtree_t;

typedef bgm::point<int, 2, bg::cs::cartesian> ipoint_t;
typedef bgm::box<ipoint_t> ibox_t;

#endif
