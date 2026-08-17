#pragma once

#include <cassert>
#include <utility>

namespace lct {

struct node {
	node* p;
	node* c[2];

	int s;

	bool flip;

	// isroot
	inline bool r() { return p == nullptr || !(this == p->c[0] || this == p->c[1]); }
	// direction
	inline bool d() { assert(!r()); return this == p->c[1]; }

	inline void update() { s = 1 + (c[0] ? c[0]->s : 0) + (c[1] ? c[1]->s : 0); }
	void propogate() {
		if(flip) {
			std::swap(c[0], c[1]);
			if(c[0]) c[0]->flip = !c[0]->flip;
			if(c[1]) c[1]->flip = !c[1]->flip;
			flip = false;
		}
	}

	// precondition: parent and current are propogated
	void rot() {
		assert(!r());

		int x = d();
		node* pa = p;
		node* ch = c[!x];

		assert(!pa->flip);
		assert(!flip);

		assert((!ch) || ch->p == this);

		if(!pa->r()) pa->p->c[pa->d()] = this;
		this->p = pa->p;

		pa->c[x] = ch;
		if(ch) ch->p = pa;

		this->c[!x] = pa;
		pa->p = this;

		pa->update();
		update();
	}

	// postcondition: always propogated
	void splay() {
		if(r()) {
			update();
			propogate();
			return;
		}

		while(!r()) {
			if(!p->r()) {
				node* gp = p->p;
				node* pa = p;
				gp->propogate();
				pa->propogate();
				propogate();
				if(d() == p->d()) {
					pa->rot();
					assert(p == pa);
				} else {
					rot();
					assert(p == gp);
				}
				rot();
			} else {
				p->propogate();
				propogate();
				rot();
				assert(r());
			}
		}
		update();
	}

	// attach on right side
	// precondition: propogated
	void make_child(node* n) {
		assert(!flip);
		assert(r());

		if(c[1]) {
			node* v = c[1];
			c[1] = nullptr;
			assert(v->r());

			update();
		}

		assert(!flip);
		assert(!c[1]);

		if(n) {

			assert(n->r());
			assert(n->p == this);

			c[1] = n;
			assert(c[1]->p == this);

			update();
		}
	}

	// postcondition: propogated
	void expose() {
		splay();
		assert(!flip);
		make_child(nullptr);
		while(p) {
			assert(r());
			p->splay();
			p->make_child(this);
			assert(!p->flip);
			assert(!flip);
			rot();
			update();
			assert(r());
		}
		assert(!p);
		assert(!c[1]);
	}

	// does not propogate
	void make_root() {
		expose();
		assert(p == nullptr);
		assert(r());
		flip = !flip;
	}

};

} // namespace lct
