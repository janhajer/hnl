#pragma once

#include <CGAL/license/Polyhedron.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Random.h>

namespace CGAL {

struct DefaultColorFunctorPolyhedron { // Default color functor; user can change it to have its own face color
    template<typename Polyhedron>
    static CGAL::Color run(const Polyhedron&, typename Polyhedron::Facet_const_handle fh) {    // use to get the mono color
        if (fh == boost::graph_traits<Polyhedron>::null_face()) {
            return CGAL::Color(100, 125, 200);       // R G B between 0-255
        }
        CGAL::Random random((unsigned int)(std::size_t)(& (*fh)));
        return get_random_color(random);
    }
};

template<class Polyhedron, class ColorFunctor>
class SimplePolyhedronViewerQt : public Basic_viewer_qt { // Viewer class for Polyhedron
    typedef Basic_viewer_qt Base;
    typedef typename Polyhedron::Traits Kernel;
    typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
    typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
public:
    SimplePolyhedronViewerQt(QWidget* parent, const char* title = "Basic Polyhedron Viewer", bool anofaces = false, const ColorFunctor& fcolor = ColorFunctor()) :
        Base(parent, title, false, true, true, true, true), // First draw: no vertex; edges, faces; mono-color; inverse normal
        m_nofaces(anofaces),
        m_fcolor(fcolor) { }
    void add(const Polyhedron& apoly) {
        polys.emplace_back(apoly);
        compute_elements();
    }
protected:
    void compute_face(Facet_const_handle fh, const Polyhedron& poly) {
        CGAL::Color c = m_fcolor.run(poly, fh);
        face_begin(c);
        Halfedge_const_handle he = fh->facet_begin();
        do {
            add_point_in_face(he->vertex()->point(), get_vertex_normal(he));
            he = he->next();
        } while (he != fh->facet_begin());
        face_end();
    }
    void compute_edge(Halfedge_const_handle he) {
        add_segment(he->vertex()->point(), he->opposite()->vertex()->point()); // We can use add_segment(p1, p2, c) with c a CGAL::Color to add a colored segment
    }
    void compute_vertex(Vertex_const_handle vh) {
        add_point(vh->point()); // We can use add_point(p, c) with c a CGAL::Color to add a colored point
    }
    void compute_elements() {
        clear();
        for (auto const& poly : polys) {
            if (!m_nofaces) for (auto f = poly.facets_begin(); f != poly.facets_end(); f++) if (f != boost::graph_traits<Polyhedron>::null_face()) compute_face(f, poly);
            for (auto e = poly.halfedges_begin(); e != poly.halfedges_end(); ++e) if (e < e->opposite()) compute_edge(e);
            for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) compute_vertex(v);
        }
    }
    virtual void keyPressEvent(QKeyEvent* e) {
        Base::keyPressEvent(e);
    }
protected:
    Local_vector get_face_normal(Halfedge_const_handle he) {
        Local_vector normal = CGAL::NULL_VECTOR;
        Halfedge_const_handle end = he;
        unsigned int nb = 0;
        do {
            internal::newell_single_step_3
            (this->get_local_point(he->vertex()->point()),
             this->get_local_point(he->next()->vertex()->point()),
             normal);
            ++nb;
            he = he->next();
        } while (he != end);
        assert(nb > 0);
        return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0 / nb));
    }
    Local_vector get_vertex_normal(Halfedge_const_handle he) {
        Local_vector normal = CGAL::NULL_VECTOR;
        Halfedge_const_handle end = he;
        do {
            if (!he->is_border()) {
                Local_vector n = get_face_normal(he);
                normal = typename Local_kernel::Construct_sum_of_vectors_3()(normal, n);
            }
            he = he->next()->opposite();
        } while (he != end);
        if (!typename Local_kernel::Equal_3()(normal, CGAL::NULL_VECTOR)) normal = (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0 / CGAL::sqrt(normal.squared_length())));
        return normal;
    }
protected:
    std::vector<Polyhedron> polys;
    bool m_nofaces;
    const ColorFunctor& m_fcolor;
};

template<class PolyhedronTraits_3, class PolyhedronItems_3, template < class T, class I, class A> class T_HDS, class Alloc>
void draw(const CGAL::Polyhedron_3 <PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>& apoly, const char* title = "Polyhedron Basic Viewer", bool nofill = false) { // Specialization of draw function.
    int argc = 1;
    const char* argv[2] = {"polyhedron_viewer", "\0"};
    QApplication app(argc, const_cast<char**>(argv));
    DefaultColorFunctorPolyhedron fcolor;
    SimplePolyhedronViewerQt<CGAL::Polyhedron_3 <PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>, DefaultColorFunctorPolyhedron>
    mainwindow(app.activeWindow(), apoly, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
}

}
