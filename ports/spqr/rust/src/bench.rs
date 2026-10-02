use spqr::spqr_tree::{PlanarSpqrTree, SpqrTree};
use spqr::spqr_tree_fast;
use spqr::spqr_tree_idiomatic as idiom;
use std::io::Read;
use std::time::Instant;

fn main() {
	let r: usize = std::env::args().nth(1).map(|s| s.parse().unwrap()).unwrap_or(5);
	let only: Option<String> = std::env::args().nth(2);
	let mut input = String::new();
	std::io::stdin().read_to_string(&mut input).unwrap();
	let mut it = input.split_ascii_whitespace().map(|x| x.parse::<i32>().unwrap());
	let nv = it.next().unwrap();
	let ne = it.next().unwrap();
	let edges: Vec<[i32; 2]> = (0..ne).map(|_| [it.next().unwrap(), it.next().unwrap()]).collect();
	let bench = |name: &str, f: &dyn Fn() -> usize| {
		if only.as_deref().is_some_and(|o| !name.starts_with(o)) {
			return;
		}
		let mut best = f64::MAX;
		let mut sink = 0;
		for _ in 0..r {
			let t0 = Instant::now();
			sink += f();
			best = best.min(t0.elapsed().as_secs_f64() * 1e3);
		}
		println!("rust {:<16} {:8.2} ms  (items={})", name, best, sink / r);
	};
	bench("spqr", &|| SpqrTree::build(nv, &edges, false, &[], &[]).par.len());
	bench("planar", &|| PlanarSpqrTree::build(nv, &edges, false, &[], &[]).par.len());
	let edges_u: Vec<[usize; 2]> = edges.iter().map(|&[u, v]| [u as usize, v as usize]).collect();
	bench("idiom spqr", &|| idiom::SpqrTree::build(nv as usize, &edges_u, false, &[], &[]).len());
	bench("idiom planar", &|| idiom::PlanarSpqrTree::build(nv as usize, &edges_u, false, &[], &[]).len());
	bench("fast spqr", &|| spqr_tree_fast::build(nv, &edges, false, &[], &[]).par.len());
	bench("fast planar", &|| spqr_tree_fast::build_planar(nv, &edges, false, &[], &[]).par.len());
}
