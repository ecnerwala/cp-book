struct PointData {
};
struct RangeData {
	RangeData(); // empty range maybe?
	explicit RangeData(const PointData&);
	friend RangeData operator + (const RangeData& a, const RangeData& b);
};

struct Op {
	Op(); // No-op
	void operator () (PointData& p);
	void operator () (RangeData& d);
	void operator () (Op& a);

	void split(RangeData& a, RangeData& b);
};

bool splitDecision(RangeData& left, RangeData& right) {
}

struct tree {
	struct node {
		node* p;
		node* c[2];
	};
};
