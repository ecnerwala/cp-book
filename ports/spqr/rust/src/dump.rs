use spqr::spqr_tree::{PlanarSpqrTree, SpqrTree};
use spqr::spqr_tree::{Csr, NodeAdj, NodeEdge, NodeType, NodeVert};
use spqr::spqr_tree_fast;
use spqr::spqr_tree_idiomatic as idiom;
use std::fmt::Write as _;
use std::io::{Read, Write};

fn main() {
	let mut input = String::new();
	std::io::stdin().read_to_string(&mut input).unwrap();
	let mut it = input.split_ascii_whitespace().map(|x| x.parse::<i32>().unwrap());
	let nv = it.next().unwrap();
	let ne = it.next().unwrap();
	let ternarize = it.next().unwrap() != 0;
	let edges: Vec<[i32; 2]> = (0..ne).map(|_| [it.next().unwrap(), it.next().unwrap()]).collect();
	let k = it.next().unwrap();
	let vert_order: Vec<i32> = (0..k).map(|_| it.next().unwrap()).collect();
	let l = it.next().unwrap();
	let edge_order: Vec<i32> = (0..l).map(|_| it.next().unwrap()).collect();

	let variant = std::env::args().nth(1).unwrap_or_else(|| "faithful".to_string());
	let (t, s) = match variant.as_str() {
		"faithful" => (
			PlanarSpqrTree::build(nv, &edges, ternarize, &vert_order, &edge_order),
			SpqrTree::build(nv, &edges, ternarize, &vert_order, &edge_order),
		),
		"fast" => (
			spqr_tree_fast::build_planar(nv, &edges, ternarize, &vert_order, &edge_order),
			spqr_tree_fast::build(nv, &edges, ternarize, &vert_order, &edge_order),
		),
		"idiomatic" => {
			let to_usize = |v: &[i32]| v.iter().map(|&x| x as usize).collect::<Vec<_>>();
			let edges_u: Vec<[usize; 2]> = edges.iter().map(|&[u, v]| [u as usize, v as usize]).collect();
			let (vo, eo) = (to_usize(&vert_order), to_usize(&edge_order));
			let p = idiom::PlanarSpqrTree::build(nv as usize, &edges_u, ternarize, &vo, &eo);
			let s = idiom::SpqrTree::build(nv as usize, &edges_u, ternarize, &vo, &eo);
			(from_idiomatic_planar(&p), from_idiomatic(&s))
		}
		_ => panic!("unknown variant {variant}"),
	};
	let mut out = String::new();
	let dump = |out: &mut String, name: &str, v: &[i32]| {
		out.push_str(name);
		out.push(':');
		for x in v {
			write!(out, " {}", x).unwrap();
		}
		out.push('\n');
	};
	dump(&mut out, "vert_index", &t.vert_index);
	dump(&mut out, "edge_index", &t.edge_index);
	dump(&mut out, "par", &t.par);
	dump(&mut out, "subtree_end", &t.subtree_end);
	out.push_str("types:");
	for x in &t.types {
		write!(out, " {}", x).unwrap();
	}
	out.push('\n');
	dump(&mut out, "orig_id", &t.orig_id);
	dump(&mut out, "ch.bounds", &t.ch.bounds);
	dump(&mut out, "ch.dat", &t.ch.dat);
	dump(&mut out, "node_verts.bounds", &t.node_verts.bounds);
	out.push_str("node_verts.dat:");
	for x in &t.node_verts.dat {
		write!(out, " {},{}", x.node, x.vert).unwrap();
	}
	out.push('\n');
	dump(&mut out, "vert_par_nv", &t.vert_par_nv);
	dump(&mut out, "node_edges.bounds", &t.node_edges.bounds);
	out.push_str("node_edges.dat:");
	for x in &t.node_edges.dat {
		write!(out, " {},{},{},{}", x.node, x.twin_ne, x.nvs[0], x.nvs[1]).unwrap();
	}
	out.push('\n');
	dump(&mut out, "node_adj.bounds", &t.node_adj.bounds);
	out.push_str("node_adj.dat:");
	for x in &t.node_adj.dat {
		write!(out, " {},{}", x.ne, x.dest_nv).unwrap();
	}
	out.push('\n');
	out.push_str("node_planar:");
	for x in &t.node_planar {
		write!(out, " {}", *x as i32).unwrap();
	}
	out.push('\n');
	dump(&mut out, "ne_rot_adj", &t.ne_rot_adj);

	writeln!(out, "nonplanar_build_same: {}", (s == t.tree) as i32).unwrap();
	std::io::stdout().write_all(out.as_bytes()).unwrap();
}

fn i(x: idiom::Idx) -> i32 {
	x.get() as i32
}
fn o(x: Option<idiom::Idx>) -> i32 {
	x.map_or(-1, i)
}
fn csr<A, B>(c: &idiom::Csr<A>, f: impl Fn(&A) -> B) -> Csr<B> {
	Csr { bounds: c.bounds.iter().map(|&b| b as i32).collect(), dat: c.dat.iter().map(f).collect() }
}
fn from_idiomatic(t: &idiom::SpqrTree) -> SpqrTree {
	let ty = |t: idiom::NodeType| match t {
		idiom::NodeType::F => NodeType::F,
		idiom::NodeType::V => NodeType::V,
		idiom::NodeType::Q => NodeType::Q,
		idiom::NodeType::I => NodeType::I,
		idiom::NodeType::O => NodeType::O,
		idiom::NodeType::S => NodeType::S,
		idiom::NodeType::P => NodeType::P,
		idiom::NodeType::R => NodeType::R,
	};
	SpqrTree {
		vert_index: t.vert_index.iter().copied().map(i).collect(),
		edge_index: t.edge_index.iter().copied().map(i).collect(),
		par: t.par.iter().copied().map(o).collect(),
		subtree_end: t.subtree_end.iter().copied().map(i).collect(),
		types: t.types.iter().copied().map(ty).collect(),
		orig_id: t.orig_id.iter().copied().map(o).collect(),
		ch: csr(&t.ch, |&x| i(x)),
		node_verts: csr(&t.node_verts, |x| NodeVert { node: i(x.node), vert: i(x.vert) }),
		vert_par_nv: t.vert_par_nv.iter().copied().map(o).collect(),
		node_edges: csr(&t.node_edges, |x| NodeEdge { node: i(x.node), twin_ne: i(x.twin_ne), nvs: x.nvs.map(i) }),
		node_adj: csr(&t.node_adj, |x| NodeAdj { ne: i(x.ne), dest_nv: i(x.dest_nv) }),
	}
}
fn from_idiomatic_planar(t: &idiom::PlanarSpqrTree) -> PlanarSpqrTree {
	PlanarSpqrTree { tree: from_idiomatic(&t.tree), node_planar: t.node_planar.clone(), ne_rot_adj: t.ne_rot_adj.iter().copied().map(o).collect() }
}
