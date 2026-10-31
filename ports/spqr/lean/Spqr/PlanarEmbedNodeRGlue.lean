import Spqr.PlanarEmbedNodeRCert

/-!
# The `R` node: the hybrid system

`skR i` is the local skeleton of the `R` node `i` relabelled to original vertices; `cornersR` lists
the corners of the executable fold (where a `V` item is hung). `rg_R` packages the skeleton, the
corner pieces and the capped children as an `RG` satisfying `RG.Hyp`, so that `RG.Hyp.glue` is the
planar witness of the node's piece.
-/

namespace Spqr.PlanarSpqrTree

open EmbedM Classical

variable (t : PlanarSpqrTree) {g : Graph} {i : Nat}

/-- The skeleton of `i` on original vertices. -/
def skR (i : Nat) : List (Nat × Nat) := mapEdges (t.sV i) (t.localSkeleton i)

theorem skR_length : (t.skR i).length = t.toSpqrTree.nEdges i := by
  simp [skR, mapEdges, localSkeleton_length]

theorem skR_getElem? (hwf : t.toSpqrTree.WF) (hi : i < t.size) {k : Nat}
    (hk : k < t.toSpqrTree.nEdges i) :
    (t.skR i)[k]? = some
      (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 - (t.toSpqrTree.nvRange i).1),
       t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 - (t.toSpqrTree.nvRange i).1)) := by
  unfold skR mapEdges
  rw [List.getElem?_map, t.localSkeleton_getElem? hwf hi hk]
  rfl

theorem vert_skR (hwf : t.toSpqrTree.WF) (hi : i < t.size) {l : Nat}
    (hl : l < 4 * t.toSpqrTree.nEdges i) :
    QE.vert (t.skR i) l = some (t.sV i
      ((if l % 4 < 2 then (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.1
        else (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2) - (t.toSpqrTree.nvRange i).1)) := by
  unfold QE.vert QE.edge QE.side
  rw [t.skR_getElem? hwf hi (k := l / 4) (by omega)]
  by_cases h : l % 4 < 2
  · rw [ite_eq_left h]; simp [show l / 2 % 2 = 0 by omega]
  · rw [ite_eq_right h]; simp [show l / 2 % 2 ≠ 0 by omega]

theorem hasEdge_sk_R (hwf : t.toSpqrTree.WF) (hrep : t.toSpqrTree.Represents g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    {k : Nat} (hk : k < t.toSpqrTree.nVerts i) : HasEdge ((t.localSkeleton i).eraseIdx 0) k := by
  have h3 : SpqrTree.ThreeConnected (t.toSpqrTree.nVerts i) (t.localSkeleton i) :=
    hrep.r_three_connected i hi hR
  have h6 := (t.shape_R hwf hi hR).2.1
  have he := t.localSkeleton_getElem? hwf hi (k := 0) (by omega)
  have hnvs := t.nvs_R hwf hi hR (k := 0) (by omega)
  rw [Nat.add_zero] at he hnvs
  exact SpqrTree.ThreeConnected.hasEdge_eraseIdx h3 hloc.verts he (by omega) hk

theorem sV_inj_R (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R) {a b : Nat} (ha : a < t.toSpqrTree.nVerts i)
    (hb : b < t.toSpqrTree.nVerts i) (h : t.sV i a = t.sV i b) : a = b := by
  have h1 := t.nvOrig_nv hwf hi ha
  have h2 := t.nvOrig_nv hwf hi hb
  unfold SpqrTree.nVerts at ha hb
  have := hsep.nv_orig_inj i ((t.toSpqrTree.nvRange i).1 + a) ((t.toSpqrTree.nvRange i).1 + b) hi
    (Or.inr (Or.inr hR)) ⟨by omega, by omega⟩ ⟨by omega, by omega⟩ (by rw [h1, h2, h])
  omega

theorem vert_lt {es : List (Nat × Nat)} {n q x : Nat} (hb : ∀ p ∈ es, p.1 < n ∧ p.2 < n)
    (h : QE.vert es q = some x) : x < n := by
  unfold QE.vert at h
  obtain ⟨p, hp, hx⟩ := Option.map_eq_some_iff.1 h
  have := hb p (List.mem_of_getElem? hp)
  split at hx <;> omega

/-- Original vertices of the node-vertices of an `R` node are vertices of `g` (each lies on a child
piece). -/
theorem sV_lt_R (hwf : t.toSpqrTree.WF) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) {k : Nat} (hk : k < t.toSpqrTree.nVerts i) :
    t.sV i k < g.nv := by
  obtain ⟨p, hp, hpk⟩ := t.hasEdge_sk_R hwf hrep hi hR hloc hk
  obtain ⟨j, hj, hj0, hpj⟩ := List.mem_eraseIdx_iff_getElem.1 hp
  rw [localSkeleton_length] at hj
  have he := t.localSkeleton_getElem? hwf hi (k := j) hj
  rw [List.getElem?_eq_getElem (by rw [localSkeleton_length]; exact hj), hpj] at he
  obtain ⟨j', rfl⟩ : ∃ j', j = j' + 1 := ⟨j - 1, by omega⟩
  obtain ⟨-, -, -, -, ρ, hc⟩ := t.child_cert_R hwf hrep.ne hsep hi hR (j := j') (by omega) s h
  obtain ⟨⟨l0, -, v0⟩, -, ⟨l2, -, v2⟩, -⟩ := hc.slots_vert
  have hb := hc.planar.verts
  cases Option.some.inj he
  rcases hpk with rfl | rfl
  · exact vert_lt hb v0
  · exact vert_lt hb v2

theorem planar_skR (hwf : t.toSpqrTree.WF) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    IsPlanarEmbedding (t.skR i) g.nv (t.nodeRot i) :=
  IsPlanarEmbedding.map hloc.verts
    (fun p hp => ⟨t.sV_lt_R hwf hrep hsep hi hR hloc s h (hloc.verts p hp).1,
      t.sV_lt_R hwf hrep hsep hi hR hloc s h (hloc.verts p hp).2⟩)
    (fun _ _ hu hv huv => t.sV_inj_R hwf hsep hi hR (hu.lt_of hloc.verts) (hv.lt_of hloc.verts) huv)
    hloc

/-- The corners of the executable fold of `i`, in order. -/
noncomputable def cornersR (i : Nat) (s : EmbedState) : List Nat :=
  (List.range (4 * t.toSpqrTree.nEdges i)).filter fun l =>
    decide (t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l)

theorem mem_cornersR {s : EmbedState} {l : Nat} :
    l ∈ t.cornersR i s ↔ t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l := by
  simp only [cornersR, List.mem_filter, List.mem_range, decide_eq_true_eq]
  exact ⟨fun h => h.2, fun h => ⟨h.2.1, h⟩⟩

theorem cornersR_nodup (s : EmbedState) : (t.cornersR i s).Nodup := List.nodup_range.filter _

theorem cornersR_corner {s : EmbedState} {v : Nat} (hv : v < (t.cornersR i s).length) :
    t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot (t.cornersR i s)[v] :=
  t.mem_cornersR.1 (List.getElem_mem hv)

theorem cornersR_inj {s : EmbedState} {v v' : Nat} (hv : v < (t.cornersR i s).length)
    (hv' : v' < (t.cornersR i s).length) (h : (t.cornersR i s)[v] = (t.cornersR i s)[v']) : v = v' :=
  (List.Nodup.getElem_inj_iff (t.cornersR_nodup s)).1 h

/-- The `v`-th corner's `V` item. -/
noncomputable def vItem (i : Nat) (s : EmbedState) (v : Nat) : Nat :=
  t.vOf (4 * (t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]!)

/-- The abstract `R`-glue data of node `i` in state `s`: skeleton on original vertices, corner `V`
pieces (with their boundary pairs located in the piece) and capped children (slots located in the
piece); the witnesses are chosen by `Classical.epsilon`. -/
noncomputable def rgR (g : Graph) (i : Nat) (s : EmbedState) : RG where
  n := g.nv
  cap := (t.skR i)[0]?.getD (0, 0)
  sk := (t.skR i).drop 1
  ρ₀ := t.nodeRot i
  nV := (t.cornersR i s).length
  ves v := (t.pieceBelow g (t.vItem i s v)).es
  vρ v := Classical.epsilon fun ρ => (t.pieceBelow g (t.vItem i s v)).OpenEmbedding s.rotAdj ρ
    (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)
    (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)
    (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 -
      (t.toSpqrTree.nvRange i).1))
  vw0 v := ((t.pieceBelow g (t.vItem i s v)).loc
    (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)).getD 0
  vw1 v := ((t.pieceBelow g (t.vItem i s v)).loc
    (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)).getD 0
  vta v := (t.cornersR i s)[v]!
  vx v := t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 -
    (t.toSpqrTree.nvRange i).1)
  ces j := (t.pieceBelow g (t.sC i)[j]!).es
  cρ j := Classical.epsilon fun ρ => (t.pieceBelow g (t.sC i)[j]!).Capped s.rotAdj ρ
    (t.rQ i s (4 * (j + 1))) (t.rQ i s (4 * (j + 1) + 1)) (t.rQ i s (4 * (j + 1) + 2))
    (t.rQ i s (4 * (j + 1) + 3))
    (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.1 - (t.toSpqrTree.nvRange i).1))
    (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.2 - (t.toSpqrTree.nvRange i).1))
  cl j k := ((t.pieceBelow g (t.sC i)[j]!).loc (t.rQ i s (4 * (j + 1) + k))).getD 0

section RGSpec

variable (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
  (g : Graph) (s : EmbedState)
include hwf hi hR

theorem rgR_m : (t.rgR g i s).m + 1 = t.toSpqrTree.nEdges i := by
  have h6 := (t.shape_R hwf hi hR).2.1
  simp only [rgR, RG.m, List.length_drop, skR_length]
  omega

theorem rgR_Sk : (t.rgR g i s).Sk = t.skR i := by
  have h6 := (t.shape_R hwf hi hR).2.1
  have hne : t.skR i ≠ [] := by
    intro h; have := t.skR_length (i := i); rw [h] at this; simp at this; omega
  obtain ⟨a, l, hl⟩ := List.exists_cons_of_ne_nil hne
  simp only [rgR, RG.Sk, hl, List.getElem?_cons_zero, Option.getD_some, List.drop_succ_cons,
    List.drop_zero]

omit hR in
theorem rgR_sk_getElem? {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) :
    (t.rgR g i s).sk[j]? = some
      (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.1 - (t.toSpqrTree.nvRange i).1),
       t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.2 - (t.toSpqrTree.nvRange i).1)) := by
  show ((t.skR i).drop 1)[j]? = _
  rw [List.getElem?_drop, Nat.add_comm, t.skR_getElem? hwf hi (by omega)]

omit hR in
theorem rgR_skE {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) :
    (t.rgR g i s).skE j =
      (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.1 - (t.toSpqrTree.nvRange i).1),
       t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.2 - (t.toSpqrTree.nvRange i).1)) := by
  unfold RG.skE; rw [t.rgR_sk_getElem? hwf hi g s hj]; rfl

theorem rgR_cap : (t.rgR g i s).cap = (t.sV i 0, t.sV i (t.toSpqrTree.nVerts i - 1)) := by
  have h6 := (t.shape_R hwf hi hR).2.1
  show (t.skR i)[0]?.getD (0, 0) = _
  rw [t.skR_getElem? hwf hi (by omega), Nat.add_zero, t.capNvs_R hwf hi hR]
  simp only [Option.getD_some, Nat.sub_self]
  rw [show (t.toSpqrTree.nvRange i).2 - 1 - (t.toSpqrTree.nvRange i).1 = t.toSpqrTree.nVerts i - 1 by
    unfold SpqrTree.nVerts; omega]

omit hwf hi hR in
theorem rgR_sep_of_mem_skR {w : Nat}
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (h : HasEdge (t.skR i) w) : ∃ k, k < t.toSpqrTree.nVerts i ∧ w = t.sV i k := by
  obtain ⟨k, hk, rfl⟩ := hasEdge_mapEdges.1 h
  exact ⟨k, hk.lt_of hloc.verts, rfl⟩

omit hi hR in
theorem hasEdge_skR_sV (hrep : t.toSpqrTree.Represents g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    {k : Nat} (hk : k < t.toSpqrTree.nVerts i) : HasEdge (t.skR i) (t.sV i k) := by
  obtain ⟨p, hp, hpk⟩ := t.hasEdge_sk_R hwf hrep hi hR hloc hk
  exact hasEdge_mapEdges.2 ⟨k, ⟨p, List.mem_of_mem_eraseIdx hp, hpk⟩, rfl⟩

variable (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g) (h : t.GluedFaces g (i + 1) s)
include hne hsep h

theorem rgR_cρ {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) :
    (t.pieceBelow g (t.sC i)[j]!).Capped s.rotAdj ((t.rgR g i s).cρ j)
      (t.rQ i s (4 * (j + 1))) (t.rQ i s (4 * (j + 1) + 1)) (t.rQ i s (4 * (j + 1) + 2))
      (t.rQ i s (4 * (j + 1) + 3)) ((t.rgR g i s).skE j).1 ((t.rgR g i s).skE j).2 := by
  rw [t.rgR_skE hwf hi g s hj]
  obtain ⟨-, -, -, -, ρ, hc⟩ := t.child_cert_R hwf hne hsep hi hR hj s h
  exact Classical.epsilon_spec (p := fun ρ => (t.pieceBelow g (t.sC i)[j]!).Capped s.rotAdj ρ
    (t.rQ i s (4 * (j + 1))) (t.rQ i s (4 * (j + 1) + 1)) (t.rQ i s (4 * (j + 1) + 2))
    (t.rQ i s (4 * (j + 1) + 3))
    (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.1 - (t.toSpqrTree.nvRange i).1))
    (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.2 - (t.toSpqrTree.nvRange i).1)))
    ⟨ρ, hc⟩

theorem rgR_cl {j k : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hk : k < 4) :
    (t.pieceBelow g (t.sC i)[j]!).loc (t.rQ i s (4 * (j + 1) + k)) = some ((t.rgR g i s).cl j k) := by
  obtain ⟨l, hl⟩ := Piece.loc_exists (t.slot_R hwf hne hsep hi hR s h hj hk).1
  show _ = some (((t.pieceBelow g (t.sC i)[j]!).loc (t.rQ i s (4 * (j + 1) + k))).getD 0)
  rw [hl]; rfl

variable (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
  (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners)
include hloc hclosed hcor

theorem rgR_vρ {v : Nat} (hv : v < (t.cornersR i s).length) :
    (t.pieceBelow g (t.vItem i s v)).OpenEmbedding s.rotAdj ((t.rgR g i s).vρ v)
      (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)
      (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) ((t.rgR g i s).vx v) := by
  have hc := t.cornersR_corner (i := i) (s := s) hv
  rw [← getElem!_pos (t.cornersR i s) v hv] at hc
  obtain ⟨-, -, a0, a1, ρ, h0, h1, ho⟩ := t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc
  have e0 : t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]! = a0 := by
    unfold PlanarSpqrTree.w0; rw [h0]; rfl
  have e1 : t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]! = a1 := by
    unfold PlanarSpqrTree.w1; rw [h1]; rfl
  exact Classical.epsilon_spec (p := fun ρ => (t.pieceBelow g (t.vItem i s v)).OpenEmbedding s.rotAdj ρ
    (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)
    (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) ((t.rgR g i s).vx v)) ⟨ρ, by rw [e0, e1]; exact ho⟩

theorem rgR_vw {v : Nat} (hv : v < (t.cornersR i s).length) :
    (t.pieceBelow g (t.vItem i s v)).loc (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) =
      some ((t.rgR g i s).vw0 v) ∧
    (t.pieceBelow g (t.vItem i s v)).loc (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) =
      some ((t.rgR g i s).vw1 v) := by
  obtain ⟨la, lb, ha, hb, -, -⟩ := (t.rgR_vρ hwf hi hR g s hne hsep h hloc hclosed hcor hv).boundary
  constructor
  · show _ = some (((t.pieceBelow g (t.vItem i s v)).loc _).getD 0); rw [ha]; rfl
  · show _ = some (((t.pieceBelow g (t.vItem i s v)).loc _).getD 0); rw [hb]; rfl

end RGSpec

section Hyp

variable (hwf : t.toSpqrTree.WF) (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
  (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
  (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
  (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) (s : EmbedState)
include hwf hi hR

omit hR in
theorem nvOf_sub {nv : Nat} (hnv : t.toSpqrTree.NvOf i nv) :
    nv - (t.toSpqrTree.nvRange i).1 < t.toSpqrTree.nVerts i ∧
    t.toSpqrTree.nvOrig nv = some (t.sV i (nv - (t.toSpqrTree.nvRange i).1)) := by
  obtain ⟨hlo, hhi⟩ := id hnv
  have hlt : nv - (t.toSpqrTree.nvRange i).1 < t.toSpqrTree.nVerts i := by
    unfold SpqrTree.nVerts; omega
  refine ⟨hlt, ?_⟩
  have := t.nvOrig_nv hwf hi hlt
  rwa [Nat.add_sub_cancel' hnv.1] at this

include hrep hloc in
theorem deg_R {l : Nat} (hl : l < 4 * t.toSpqrTree.nEdges i) : (t.nodeRot i).rot l / 4 ≠ l / 4 := by
  have hs : (t.nodeRot i).size = 4 * t.toSpqrTree.nEdges i := by
    rw [hloc.size, t.localSkeleton_length]
  have he := t.localSkeleton_getElem? hwf hi (k := l / 4) (by omega)
  have hnvs := t.nvs_R hwf hi hR (k := l / 4) (by omega)
  have huv : (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.1 - (t.toSpqrTree.nvRange i).1 ≠
      (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 - (t.toSpqrTree.nvRange i).1 := by omega
  have h3 : SpqrTree.ThreeConnected (t.toSpqrTree.nVerts i) (t.localSkeleton i) :=
    hrep.r_three_connected i hi hR
  have hconn := SpqrTree.ThreeConnected.edgesConn_eraseIdx h3 hloc.verts he huv
  have hget : (t.nodeRot i).get (4 * (l / 4) + l % 4) = some ((t.nodeRot i).rot l) := by
    rw [show 4 * (l / 4) + l % 4 = l by omega]
    exact RotationSystem.get_eq_rot hloc.total (by rw [hs]; exact hl)
  exact deg_of_edgesConn hloc he huv hconn (l % 4) (by omega) _ hget

include hloc hclosed hcor in
/-- The `v`-th corner: its quarter-edge, its node-vertex (inner) and its `V` item. -/
theorem corner_data {v : Nat} (hv : v < (t.cornersR i s).length) :
    t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot (t.cornersR i s)[v]! ∧
    t.toSpqrTree.NvOf i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 ∧
    ¬ t.toSpqrTree.CapEnd i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 ∧
    t.vItem i s v = t.sW (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 ∧
    t.CornerAt i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2
      (4 * (t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]!) := by
  have hc := t.cornersR_corner (i := i) (s := s) hv
  rw [← getElem!_pos (t.cornersR i s) v hv] at hc
  obtain ⟨hnv, hnce, -⟩ := t.corner_nv_R hwf hi hR hloc hclosed hcor s hc
  exact ⟨hc, hnv, hnce, (t.vOf_R hwf hi hR hc.2.1).1, t.corner_R hwf hi hR hloc hclosed s hc⟩

include hloc hclosed hcor in
theorem corner_nv_inj {v v' : Nat} (hv : v < (t.cornersR i s).length)
    (hv' : v' < (t.cornersR i s).length)
    (heq : (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 =
      (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v']! / 4]!).nvs.2) : v = v' := by
  obtain ⟨-, hnv, hnce, -, hca⟩ := t.corner_data hwf hi hR hloc hclosed hcor s hv
  obtain ⟨-, -, -, -, hca'⟩ := t.corner_data hwf hi hR hloc hclosed hcor s hv'
  rw [← heq] at hca'
  obtain ⟨ta, -, -, huniq⟩ := hcor.inner i _ hi (Or.inr (Or.inr hR)) hnv hnce
  have e1 := huniq _ hca
  have e2 := huniq _ hca'
  apply t.cornersR_inj hv hv'
  rw [← getElem!_pos (t.cornersR i s) v hv, ← getElem!_pos (t.cornersR i s) v' hv']
  omega

include hsep hloc hclosed hcor in
theorem vItem_inj {v v' : Nat} (hv : v < (t.cornersR i s).length)
    (hv' : v' < (t.cornersR i s).length) (heq : t.vItem i s v = t.vItem i s v') : v = v' := by
  obtain ⟨-, hnv, hnce, hw, -⟩ := t.corner_data hwf hi hR hloc hclosed hcor s hv
  obtain ⟨-, hnv', hnce', hw', -⟩ := t.corner_data hwf hi hR hloc hclosed hcor s hv'
  obtain ⟨⟨d, hd, hdv⟩, -, -, -⟩ := t.sW_R (g := g) hwf hsep hi hR hnv hnce
  obtain ⟨⟨d', hd', hdv'⟩, -, -, -⟩ := t.sW_R (g := g) hwf hsep hi hR hnv' hnce'
  refine t.corner_nv_inj hwf hi hR hloc hclosed hcor s hv hv' ?_
  exact hsep.nv_vert_inj i _ _ d d' hi (Or.inr (Or.inr hR)) hnv hnv' hd hd'
    (by rw [hdv, hdv', ← hw, ← hw', heq])

include hsep hloc hclosed hcor in
/-- A corner's `V` item is attached only at its node-vertex. -/
theorem nvInc_V_R {v nv : Nat} (hv : v < (t.cornersR i s).length) (hnv : t.toSpqrTree.NvOf i nv)
    (hinc : t.toSpqrTree.NvInc i (t.vItem i s v) nv) :
    nv = (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 := by
  obtain ⟨-, hnv', hnce', hw, -⟩ := t.corner_data hwf hi hR hloc hclosed hcor s hv
  obtain ⟨⟨d', hd', hdv'⟩, -, htV, -⟩ := t.sW_R (g := g) hwf hsep hi hR hnv' hnce'
  rcases hinc with ⟨d, hd, hdv⟩ | ⟨ne, tw, d, -, -, hcne, -, -⟩
  · exact hsep.nv_vert_inj i nv _ d d' hi (Or.inr (Or.inr hR)) hnv hnv' hd hd'
      (by rw [hdv, hdv', hw])
  · exfalso
    unfold SpqrTree.capNe SpqrTree.hasCap at hcne
    rw [hw, htV] at hcne
    simp [NodeType.isNode] at hcne

include hsep in
/-- A capped child of an `R` node is attached only at the endpoints of its node-edge. -/
theorem nvInc_C_R {j nv : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hinc : t.toSpqrTree.NvInc i (t.sC i)[j]! nv) :
    nv = (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.1 ∨
    nv = (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.2 := by
  have hc := t.sC_getElem?_R hwf hi hR hj
  obtain ⟨-, hcV, hcap, -, -, htw, -, -⟩ := t.child_R (g := g) hwf hsep hi hR hj hc
  rcases hinc with ⟨d, hd, hdv⟩ | ⟨ne, tw, d, -, htw', hcne, hd, hend⟩
  · exfalso
    have := hwf.own.nv_vert nv (Array.getElem?_eq_some_iff.1 hd).1
    rw [getElem!_of_getElem? hd, hdv] at this
    exact hcV this
  · have hcne' : t.toSpqrTree.capNe (t.sC i)[j]! = some (t.toSpqrTree.neRange (t.sC i)[j]!).1 := by
      simp only [SpqrTree.capNe, hcap, ↓reduceIte]
    rw [hcne'] at hcne
    cases hcne
    have h1 := hwf.twins.twin_invol _ _ htw'
    have h2 := hwf.twins.twin_invol _ _ htw
    rw [h1] at h2
    cases h2
    rw [getElem!_of_getElem? hd]
    omega

include hrep hsep hloc in
/-- Two distinct children of the `R` node meet only on skeleton vertices. -/
theorem attach_skR {a b w : Nat} (ha : a ∈ t.children i) (hb : b ∈ t.children i) (hab : a ≠ b)
    (hta : t.toSpqrTree.Touches g a w) (htb : t.toSpqrTree.Touches g b w) : HasEdge (t.skR i) w := by
  obtain ⟨nv, hnv, horig⟩ := hsep.node_attach i a b w hi (Or.inr (Or.inr hR)) (by rw [t.children_eq]; exact ha)
    (by rw [t.children_eq]; exact hb) hab hta htb
  obtain ⟨hlt, horig'⟩ := t.nvOf_sub hwf hi hnv
  rw [horig'] at horig
  cases horig
  exact t.hasEdge_skR_sV hwf g hrep hi hR hloc hlt

omit hwf hi hR in
theorem lt_size_of_get {ρ : RotationSystem} {q r : Nat} (h : ρ.get q = some r) : q < ρ.size := by
  obtain ⟨a, ha, -⟩ := Option.bind_eq_some_iff.1 h
  exact (Array.getElem?_eq_some_iff.1 ha).1

include hrep hsep hloc hclosed hcor in
theorem hyp_R (hne : t.ne = g.ne) (h : t.GluedFaces g (i + 1) s) : (t.rgR g i s).Hyp := by
  have hSk := t.rgR_Sk hwf hi hR g s
  have hm := t.rgR_m hwf hi hR g s
  have hcap := t.rgR_cap hwf hi hR g s
  have h4 := (hrep.r_three_connected i hi hR).1
  have h6 := (t.shape_R hwf hi hR).2.1
  have hC := fun {j} (hj : j < (t.rgR g i s).m) =>
    t.rgR_cρ hwf hi hR g s hne hsep h (j := j) (by omega)
  have hcl := fun {j k} (hj : j < (t.rgR g i s).m) (hk : k < 4) =>
    t.rgR_cl hwf hi hR g s hne hsep h (j := j) (k := k) (by omega) hk
  have hV := fun {v} (hv : v < (t.rgR g i s).nV) =>
    t.rgR_vρ hwf hi hR g s hne hsep h hloc hclosed hcor (v := v) hv
  have hvw := fun {v} (hv : v < (t.rgR g i s).nV) =>
    t.rgR_vw hwf hi hR g s hne hsep h hloc hclosed hcor (v := v) hv
  have hcd := fun {v} (hv : v < (t.rgR g i s).nV) =>
    t.corner_data hwf hi hR hloc hclosed hcor s (v := v) hv
  have hvget : ∀ v, v < (t.rgR g i s).nV →
      ((t.rgR g i s).vρ v).get ((t.rgR g i s).vw0 v) = some ((t.rgR g i s).vw1 v) ∧
      QE.vert ((t.rgR g i s).ves v) ((t.rgR g i s).vw0 v) = some ((t.rgR g i s).vx v) := by
    intro v hv
    obtain ⟨la, lb, ha, hb, hget, hvert⟩ := (hV hv).boundary
    obtain ⟨e0, e1⟩ := hvw hv
    rw [e0] at ha; rw [e1] at hb
    cases ha; cases hb
    exact ⟨hget, hvert⟩
  have hcget : ∀ j, j < (t.rgR g i s).m →
      ((t.rgR g i s).cρ j).get ((t.rgR g i s).cl j 0) = some ((t.rgR g i s).cl j 1) ∧
      ((t.rgR g i s).cρ j).get ((t.rgR g i s).cl j 2) = some ((t.rgR g i s).cl j 3) ∧
      QE.vert ((t.rgR g i s).ces j) ((t.rgR g i s).cl j 0) = some ((t.rgR g i s).skE j).1 ∧
      QE.vert ((t.rgR g i s).ces j) ((t.rgR g i s).cl j 2) = some ((t.rgR g i s).skE j).2 ∧
      ((t.rgR g i s).cρ j).SameFaceOrbit ((t.rgR g i s).cl j 0) ((t.rgR g i s).cl j 2) := by
    intro j hj
    obtain ⟨l0, l1, h0, h1, g0, v0⟩ := (hC hj).pair0
    obtain ⟨l2, l3, h2, h3, g2, v2⟩ := (hC hj).pair2
    obtain ⟨l0', l2', h0', h2', hf⟩ := (hC hj).face
    have e0 := hcl hj (by omega : 0 < 4)
    rw [Nat.add_zero] at e0
    have e1 := hcl hj (by omega : 1 < 4)
    have e2 := hcl hj (by omega : 2 < 4)
    have e3 := hcl hj (by omega : 3 < 4)
    rw [e0] at h0 h0'; rw [e1] at h1; rw [e2] at h2 h2'; rw [e3] at h3
    cases h0; cases h1; cases h2; cases h3; cases h0'; cases h2'
    exact ⟨g0, g2, v0, v2, hf⟩
  have hsk : ∀ j p, (t.rgR g i s).sk[j]? = some p → j < (t.rgR g i s).m ∧ p = (t.rgR g i s).skE j := by
    intro j p hp
    exact ⟨(List.getElem?_eq_some_iff.1 hp).1, by unfold RG.skE; rw [hp]; rfl⟩
  have hmem_skR : ∀ w, HasEdge (t.rgR g i s).Sk w → ∃ k, k < t.toSpqrTree.nVerts i ∧ w = t.sV i k := by
    intro w hw
    rw [hSk] at hw
    exact t.rgR_sep_of_mem_skR hloc hw
  have hcm : ∀ j, j < (t.rgR g i s).m → (t.sC i)[j]! ∈ t.children i ∧ t.toSpqrTree.type (t.sC i)[j]! ≠ .V :=
    fun j hj => ⟨(t.child_cert_R hwf hne hsep hi hR (by omega) s h).1,
      (t.child_R (g := g) hwf hsep hi hR (by omega) (t.sC_getElem?_R hwf hi hR (by omega))).2.1⟩
  have hvm : ∀ v, v < (t.rgR g i s).nV → t.vItem i s v ∈ t.children i ∧ t.toSpqrTree.type (t.vItem i s v) = .V := by
    intro v hv
    obtain ⟨-, hnv, hnce, hw, -⟩ := hcd hv
    obtain ⟨-, hmem, htV, -⟩ := t.sW_R (g := g) hwf hsep hi hR hnv hnce
    rw [hw]; exact ⟨hmem, htV⟩
  have hsub : ∀ nv, t.toSpqrTree.NvOf i nv → nv - (t.toSpqrTree.nvRange i).1 < t.toSpqrTree.nVerts i :=
    fun nv hnv => (t.nvOf_sub hwf hi hnv).1
  have hinj := fun {a b} (ha : a < t.toSpqrTree.nVerts i) (hb : b < t.toSpqrTree.nVerts i) =>
    t.sV_inj_R hwf hsep hi hR ha hb
  refine
    { planar₀ := by rw [hSk]; exact t.planar_skR hwf hrep hsep hi hR hloc s h
      v_planar := fun v hv => (hV hv).planar
      v_get := fun v hv => (hvget v hv).1
      v_lt := fun v hv => lt_size_of_get (hvget v hv).1
      v_even := fun v hv => by
        obtain ⟨e0, -⟩ := hvw hv
        rw [Piece.loc_mod_two e0]; exact (hV hv).left_dir
      ta_ge := fun v hv => (hcd hv).1.1
      ta_lt := fun v hv => by
        show (t.cornersR i s)[v]! < 4 * ((t.rgR g i s).m + 1)
        have := (hcd hv).1.2.1; omega
      ta_even := fun v hv => by
        show (t.cornersR i s)[v]! % 2 = 0
        have := (hcd hv).1.2.2.1; omega
      ta_vert := fun v hv => by
        have hc := (hcd hv).1
        show QE.vert (t.rgR g i s).Sk (t.cornersR i s)[v]! = _
        rw [hSk, t.vert_skR hwf hi hc.2.1, ite_eq_right (by have := hc.2.2.1; omega)]
        rfl
      w_vert := fun v hv => (hvget v hv).2
      x_inj := fun v v' hv hv' heq => by
        obtain ⟨-, hnv, -, -, -⟩ := hcd hv
        obtain ⟨-, hnv', -, -, -⟩ := hcd hv'
        have heq' : t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v]! / 4]!).nvs.2 -
            (t.toSpqrTree.nvRange i).1) =
          t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (t.cornersR i s)[v']! / 4]!).nvs.2 -
            (t.toSpqrTree.nvRange i).1) := heq
        have := hinj (hsub _ hnv) (hsub _ hnv') heq'
        obtain ⟨hl1, hl2⟩ := id hnv; obtain ⟨hl1', hl2'⟩ := id hnv'
        exact t.corner_nv_inj hwf hi hR hloc hclosed hcor s hv hv' (by omega)
      x_cap := fun v hv => by
        obtain ⟨-, hnv, hnce, -, -⟩ := hcd hv
        rw [hcap]
        obtain ⟨hlo, hhi⟩ := id hnv
        constructor
        · intro heq
          have := hinj (hsub _ hnv) (by omega) heq
          exact hnce ((t.capEnd_R hwf hi hR _).2 (Or.inl (by omega)))
        · intro heq
          have := hinj (hsub _ hnv) (by omega) heq
          unfold SpqrTree.nVerts at this
          exact hnce ((t.capEnd_R hwf hi hR _).2 (Or.inr (by omega)))
      v_sep := fun v w hv hw hwv => by
        obtain ⟨k, hk, rfl⟩ := hmem_skR w hw
        obtain ⟨-, hnv, -, -, -⟩ := hcd hv
        have hnv' : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + k) :=
          ⟨by omega, by unfold SpqrTree.nVerts at hk; omega⟩
        have hinc := hsep.node_touch i _ _ _ hi (Or.inr (Or.inr hR))
          (by rw [t.children_eq]; exact (hvm v hv).1)
          (t.touches_of_hasEdge hwf g hne hwv) hnv' (t.nvOrig_nv hwf hi hk)
        have := t.nvInc_V_R hwf hsep hi hR hloc hclosed hcor s hv hnv' hinc
        show t.sV i k = t.sV i (_ - (t.toSpqrTree.nvRange i).1)
        rw [← this, Nat.add_sub_cancel_left]
      vv_sep := fun v v' w hv hv' hne' hw hw' => by
        rw [hSk]
        exact t.attach_skR hwf hrep hsep hi hR hloc (hvm v hv).1 (hvm v' hv').1
          (fun heq => hne' (t.vItem_inj hwf hsep hi hR hloc hclosed hcor s hv hv' heq))
          (t.touches_of_hasEdge hwf g hne hw) (t.touches_of_hasEdge hwf g hne hw')
      c_planar := fun j hj => (hC hj).planar
      c_get0 := fun j hj => (hcget j hj).1
      c_get2 := fun j hj => (hcget j hj).2.1
      c_even0 := fun j hj => by
        have e0 := hcl hj (by omega : 0 < 4)
        rw [Nat.add_zero] at e0
        rw [Piece.loc_mod_two e0]; exact (hC hj).dir0
      c_even2 := fun j hj => by
        rw [Piece.loc_mod_two (hcl hj (by omega : 2 < 4))]; have := (hC hj).dir2; omega
      c_lt0 := fun j hj => lt_size_of_get (hcget j hj).1
      c_lt2 := fun j hj => lt_size_of_get (hcget j hj).2.1
      c_vert := fun j p hp => by
        obtain ⟨hj, rfl⟩ := hsk j p hp
        refine ⟨(hcget j hj).2.2.1, (hcget j hj).2.2.2.1, ?_⟩
        rw [t.rgR_skE hwf hi g s (by omega)]
        exact t.sV_ne_R hwf hsep hi hR (k := j + 1) (by omega)
      c_face := fun j hj => by
        have hpl := (hC hj).planar
        exact (RotationSystem.sameFaceOrbit_iff hpl.total hpl.involution hpl.size
          (lt_size_of_get (hcget j hj).1)).1 (hcget j hj).2.2.2.2
      c_sep := fun j p w hp hwc hwS => by
        obtain ⟨hj, rfl⟩ := hsk j p hp
        obtain ⟨k, hk, rfl⟩ := hmem_skR w hwS
        have hnv' : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + k) :=
          ⟨by omega, by unfold SpqrTree.nVerts at hk; omega⟩
        have hinc := hsep.node_touch i _ _ _ hi (Or.inr (Or.inr hR))
          (by rw [t.children_eq]; exact (hcm j hj).1)
          (t.touches_of_hasEdge hwf g hne hwc) hnv' (t.nvOrig_nv hwf hi hk)
        rw [t.rgR_skE hwf hi g s (by omega)]
        rcases t.nvInc_C_R hwf hsep hi hR (by omega) hinc with e | e
        · left; show t.sV i k = t.sV i (_ - (t.toSpqrTree.nvRange i).1); rw [← e, Nat.add_sub_cancel_left]
        · right; show t.sV i k = t.sV i (_ - (t.toSpqrTree.nvRange i).1); rw [← e, Nat.add_sub_cancel_left]
      cc_sep := fun j j' w hj hj' hjj hw hw' => by
        rw [hSk]
        exact t.attach_skR hwf hrep hsep hi hR hloc (hcm j hj).1 (hcm j' hj').1
          (fun heq => hjj (t.sC_inj_R hwf hi hR (by omega) (by omega) heq))
          (t.touches_of_hasEdge hwf g hne hw) (t.touches_of_hasEdge hwf g hne hw')
      cv_sep := fun j v w hj hv hw hw' => by
        rw [hSk]
        exact t.attach_skR hwf hrep hsep hi hR hloc (hcm j hj).1 (hvm v hv).1
          (fun heq => (hcm j hj).2 (heq ▸ (hvm v hv).2))
          (t.touches_of_hasEdge hwf g hne hw) (t.touches_of_hasEdge hwf g hne hw')
      deg₀ := fun l hl => t.deg_R hwf hrep hi hR hloc (by omega)
      cap_ne := by
        rw [hcap]
        intro heq
        have := hinj (by omega) (by omega) heq
        omega
      cap_conn := ?_ }
  rw [hcap]
  show EdgesConn ((t.skR i).drop 1) _ _
  have hdrop : (t.skR i).drop 1 = mapEdges (t.sV i) ((t.localSkeleton i).eraseIdx 0) := by
    simp [skR, mapEdges, List.drop_one, List.eraseIdx_zero, List.map_tail]
  rw [hdrop]
  have hb : ∀ p ∈ (t.localSkeleton i).eraseIdx 0, p.1 < t.toSpqrTree.nVerts i ∧ p.2 < t.toSpqrTree.nVerts i :=
    fun p hp => hloc.verts p (List.mem_of_mem_eraseIdx hp)
  have hinj' : ∀ u v, HasEdge ((t.localSkeleton i).eraseIdx 0) u → HasEdge ((t.localSkeleton i).eraseIdx 0) v →
      t.sV i u = t.sV i v → u = v :=
    fun _ _ hu hv huv => hinj (hu.lt_of hb) (hv.lt_of hb) huv
  have hx := t.hasEdge_sk_R hwf hrep hi hR hloc (k := 0) (by omega)
  have h3 : SpqrTree.ThreeConnected (t.toSpqrTree.nVerts i) (t.localSkeleton i) :=
    hrep.r_three_connected i hi hR
  have he := t.localSkeleton_getElem? hwf hi (k := 0) (by omega)
  rw [Nat.add_zero, t.capNvs_R hwf hi hR] at he
  simp only [Nat.sub_self] at he
  have hconn := SpqrTree.ThreeConnected.edgesConn_eraseIdx h3 hloc.verts he (by
    have := t.nvs_R hwf hi hR (k := 0) (by omega)
    rw [Nat.add_zero, t.capNvs_R hwf hi hR] at this
    omega)
  rw [show t.toSpqrTree.nVerts i - 1 = (t.toSpqrTree.nvRange i).2 - 1 - (t.toSpqrTree.nvRange i).1 by
    unfold SpqrTree.nVerts; omega]
  exact (edgesConn_mapEdges_iff (f := t.sV i) hinj' hx).2 ⟨_, hconn, rfl⟩

end Hyp

end Spqr.PlanarSpqrTree
