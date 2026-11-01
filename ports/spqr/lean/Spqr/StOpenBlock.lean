import Spqr.StSim
import Spqr.StSimLemmas
/-!
# The open block of the truncated reference

`openBlock g fs dirs ps` is the block that the stack segment of the subtree at the bottom of the
open path `fs` belongs to, when that subtree has produced the pieces `ps`: the pieces are plugged
into the frames of `fs` from the bottom up (`retPieces`, the returning-edge assembly of `refOut`)
until the deepest block-boundary frame, whose edge is the root of the block (`none` at a DFS root).
It is the last block of `refBlocks g [truncTree fs t]` for the subtree `t` with pieces `ps`.
-/

namespace Spqr

/-- The pieces of a returning out-edge `o` (`o.cls.lowval d < d`) of `v`, given the pieces
`childPs` of its child subtree (`refOut`'s returning branch). -/
def retPieces (g : Graph) (v d : Nat) (dirs : List Bool) (o : DfsOut) (hv : Bool)
    (childPs : List StPiece) : List StPiece :=
  let lowDir := dirs.getD (o.cls.lowval d) false
  let sd := !lowDir
  let pre := if !hv && o.cls.isType1 then [StPiece.mk sd [vertItem v]] else []
  let hv₁ := (!hv && o.cls.isType1) || hv
  let mid := match o with
    | .tree e _ _ =>
      let sub := childPs ++ [StPiece.mk sd [edgeItem g e]]
      if hv₁ then [StPiece.mk lowDir (stNest sub)] else sub
    | .back e _ _ => [StPiece.mk lowDir [edgeItem g e]]
  let post := if hv₁ then [] else [StPiece.mk sd [vertItem v]]
  pre ++ mid ++ post

/-- The pieces of the child subtree of `o` (`[]` for a back edge). -/
def childPieces (g : Graph) (d : Nat) (dirs : List Bool) (o : DfsOut) : List StPiece :=
  match o with
  | .tree _ _ child => (refTree g child (d + 1) (dirs ++ [!dirs.getD (o.cls.lowval d) false])).1
  | .back .. => []

theorem refOut_ret_pieces {g : Graph} {v d : Nat} {dirs : List Bool} {o : DfsOut} {hv : Bool}
    (h : ¬ d ≤ o.cls.lowval d) :
    (refOut g v d dirs o hv).1 = retPieces g v d dirs o hv (childPieces g d dirs o) := by
  cases o with
  | tree e cls child =>
    have h' : ¬ d ≤ cls.lowval d := h
    cases hv <;> cases hT : cls.isType1 <;>
      simp [refOut, retPieces, childPieces, DfsOut.cls, h', hT]
  | back e dest cls =>
    have h' : ¬ d ≤ cls.lowval d := h
    cases hv <;> cases hT : cls.isType1 <;>
      simp [refOut, retPieces, childPieces, DfsOut.cls, h', hT]

/-- Bottom-up plugging of the pieces `ps` into the frames (frame `k` at depth `k`). -/
def openBlockP (g : Graph) (dirs : List Bool) :
    Nat → List PathFrame → List StPiece → Option (Nat × Nat) × List StPiece
  | _, [], ps => (none, ps)
  | k, f :: fs, ps =>
    match openBlockP g dirs (k + 1) fs ps with
    | (some r, ps') => (some r, ps')
    | (none, ps') =>
      if k ≤ f.o.cls.lowval k then (some (f.v, f.o.dest), ps')
      else (none, (refOuts g f.v k (dirs.take k) f.done false).1 ++
        retPieces g f.v k (dirs.take k) f.o (refOuts g f.v k (dirs.take k) f.done false).2.2 ps')

/-- The open block under the path `fs` (`dirs` the directions along it) for bottom pieces `ps`. -/
def openBlock (g : Graph) (fs : List PathFrame) (dirs : List Bool) (ps : List StPiece) : StBlock :=
  ⟨(openBlockP g dirs 0 fs ps).1, stNest (openBlockP g dirs 0 fs ps).2⟩

theorem openBlockP_nil (g : Graph) (dirs : List Bool) (k : Nat) (ps : List StPiece) :
    openBlockP g dirs k [] ps = (none, ps) := rfl

theorem openBlockP_snoc_bd (g : Graph) (dirs : List Bool) (f : PathFrame) (ps : List StPiece) :
    ∀ (fs : List PathFrame) (k : Nat), k + fs.length ≤ f.o.cls.lowval (k + fs.length) →
    openBlockP g dirs k (fs ++ [f]) ps = (some (f.v, f.o.dest), ps)
  | [], k, h => by
    simp only [List.nil_append, List.length_nil, Nat.add_zero] at h ⊢
    simp [openBlockP, h]
  | f' :: fs, k, h => by
    rw [List.cons_append, openBlockP, openBlockP_snoc_bd g dirs f ps fs (k + 1) (by
      simp only [List.length_cons] at h; rwa [Nat.add_assoc, Nat.add_comm 1])]

theorem openBlockP_snoc_ret (g : Graph) (dirs : List Bool) (f : PathFrame) (ps : List StPiece) :
    ∀ (fs : List PathFrame) (k : Nat), ¬ k + fs.length ≤ f.o.cls.lowval (k + fs.length) →
    openBlockP g dirs k (fs ++ [f]) ps = openBlockP g dirs k fs
      ((refOuts g f.v (k + fs.length) (dirs.take (k + fs.length)) f.done false).1 ++
        retPieces g f.v (k + fs.length) (dirs.take (k + fs.length)) f.o
          (refOuts g f.v (k + fs.length) (dirs.take (k + fs.length)) f.done false).2.2 ps)
  | [], k, h => by
    simp only [List.nil_append, List.length_nil, Nat.add_zero] at h ⊢
    simp [openBlockP, h]
  | f' :: fs, k, h => by
    have e : k + (f' :: fs).length = (k + 1) + fs.length := by simp only [List.length_cons]; omega
    rw [e] at h ⊢
    rw [List.cons_append, openBlockP, openBlockP_snoc_ret g dirs f ps fs (k + 1) h, openBlockP]

theorem openBlockP_congr (g : Graph) {dirs dirs' : List Bool} :
    ∀ (fs : List PathFrame) (k : Nat) (ps : List StPiece),
    (∀ j, j < k + fs.length → dirs.take j = dirs'.take j) →
    openBlockP g dirs k fs ps = openBlockP g dirs' k fs ps
  | [], _, _, _ => rfl
  | f :: fs, k, ps, h => by
    rw [openBlockP, openBlockP, openBlockP_congr g fs (k + 1) ps (fun j hj => h j (by
      simp only [List.length_cons]; omega)), h k (by simp only [List.length_cons]; omega)]

theorem openBlock_snoc_bd {g : Graph} {fs : List PathFrame} {dirs : List Bool} {f : PathFrame}
    {ps : List StPiece} (h : fs.length ≤ f.o.cls.lowval fs.length) :
    openBlock g (fs ++ [f]) dirs ps = ⟨some (f.v, f.o.dest), stNest ps⟩ := by
  unfold openBlock
  rw [openBlockP_snoc_bd g dirs f ps fs 0 (by simpa using h)]

theorem openBlock_snoc_ret {g : Graph} {fs : List PathFrame} {dirs : List Bool} {f : PathFrame}
    {ps : List StPiece} {x : Bool} (hd : dirs.length = fs.length)
    (h : ¬ fs.length ≤ f.o.cls.lowval fs.length) :
    openBlock g (fs ++ [f]) (dirs ++ [x]) ps = openBlock g fs dirs
      ((refOuts g f.v fs.length dirs f.done false).1 ++
        retPieces g f.v fs.length dirs f.o (refOuts g f.v fs.length dirs f.done false).2.2 ps) := by
  unfold openBlock
  have h1 : (dirs ++ [x]).take fs.length = dirs := by
    rw [← hd, List.take_append_of_le_length (Nat.le_refl _), List.take_length]
  have h2 : ∀ j, j < 0 + fs.length → (dirs ++ [x]).take j = dirs.take j := fun j hj =>
    List.take_append_of_le_length (by omega)
  rw [openBlockP_snoc_ret g _ f ps fs 0 (by simpa using h), Nat.zero_add, h1,
    openBlockP_congr g fs 0 _ h2]

theorem stNest_append' (ps qs : List StPiece) :
    stNest (ps ++ qs) = stNestL qs ++ stNest ps ++ stNestR qs := by
  simp [stNest, stNestL_append, stNestR_append]

theorem stNestL_single (b : Bool) (L : List ItemId) : stNestL [⟨b, L⟩] = if b then [] else L := by
  simp [stNestL]

theorem stNestR_single (b : Bool) (L : List ItemId) : stNestR [⟨b, L⟩] = if b then L else [] := by
  simp [stNestR]

/-- The st-order of a returning tree edge's frame is a context around the child's st-order. -/
theorem frame_ctx (g : Graph) (v d : Nat) (dirs : List Bool) (e : Nat) (cls : OutClass)
    (child : DfsTree) (hv : Bool) (qs : List StPiece) (hq : hv = false → qs = []) :
    ∃ A B, ∀ ps : List StPiece,
      stNest (qs ++ retPieces g v d dirs (.tree e cls child) hv ps) = A ++ stNest ps ++ B := by
  generalize hl : dirs.getD (cls.lowval d) false = ℓ
  cases hv with
  | true =>
    cases ℓ with
    | false =>
      refine ⟨[], [edgeItem g e] ++ stNest qs, fun ps => ?_⟩
      simp only [retPieces, DfsOut.cls, hl]
      simp [stNest_append', stNestL_single, stNestR_single]
    | true =>
      refine ⟨stNest qs ++ [edgeItem g e], [], fun ps => ?_⟩
      simp only [retPieces, DfsOut.cls, hl]
      simp [stNest_append', stNestL_single, stNestR_single]
  | false =>
    have hq' := hq rfl
    subst hq'
    cases hT : cls.isType1 with
    | true =>
      cases ℓ with
      | false =>
        refine ⟨[], [edgeItem g e, vertItem v], fun ps => ?_⟩
        simp only [retPieces, DfsOut.cls, hl, hT]
        simp [stNest, stNestL_append, stNestR_append, stNestL, stNestR]
      | true =>
        refine ⟨[vertItem v, edgeItem g e], [], fun ps => ?_⟩
        simp only [retPieces, DfsOut.cls, hl, hT]
        simp [stNest, stNestL_append, stNestR_append, stNestL, stNestR]
    | false =>
      cases ℓ with
      | false =>
        refine ⟨[], [edgeItem g e, vertItem v], fun ps => ?_⟩
        simp only [retPieces, DfsOut.cls, hl, hT]
        simp [stNest, stNestL_append, stNestR_append, stNestL, stNestR]
      | true =>
        refine ⟨[vertItem v, edgeItem g e], [], fun ps => ?_⟩
        simp only [retPieces, DfsOut.cls, hl, hT]
        simp [stNest, stNestL_append, stNestR_append, stNestL, stNestR]

/-- Frames whose out-edges are tree edges. -/
def TreeFrames (fs : List PathFrame) : Prop := ∀ f ∈ fs, ∃ e cls child, f.o = .tree e cls child

theorem openBlockP_ctx (g : Graph) (dirs : List Bool) :
    ∀ (fs : List PathFrame) (k : Nat), TreeFrames fs →
    ∃ r A B, ∀ ps, (openBlockP g dirs k fs ps).1 = r ∧
      stNest (openBlockP g dirs k fs ps).2 = A ++ stNest ps ++ B
  | [], _, _ => ⟨none, [], [], fun ps => ⟨rfl, by simp [openBlockP]⟩⟩
  | f :: fs, k, hT => by
    obtain ⟨r', A', B', ih⟩ := openBlockP_ctx g dirs fs (k + 1) fun f' hf' => hT f' (List.mem_cons_of_mem _ hf')
    obtain ⟨e, cls, child, ho⟩ := hT f List.mem_cons_self
    cases r' with
    | some r' =>
      refine ⟨some r', A', B', fun ps => ?_⟩
      obtain ⟨h1, h2⟩ := ih ps
      rcases hx : openBlockP g dirs (k + 1) fs ps with ⟨r₀, ps₀⟩
      rw [hx] at h1 h2
      simp only at h1 h2
      subst h1
      simp only [openBlockP, hx]
      exact ⟨by first | trivial | rfl, h2⟩
    | none =>
      by_cases hk : k ≤ f.o.cls.lowval k
      · refine ⟨some (f.v, f.o.dest), A', B', fun ps => ?_⟩
        obtain ⟨h1, h2⟩ := ih ps
        rcases hx : openBlockP g dirs (k + 1) fs ps with ⟨r₀, ps₀⟩
        rw [hx] at h1 h2
        simp only at h1 h2
        subst h1
        simp only [openBlockP, hx, hk, ↓reduceIte]
        exact ⟨by first | trivial | rfl, h2⟩
      · obtain ⟨A'', B'', hf⟩ := frame_ctx g f.v k (dirs.take k) e cls child
          (refOuts g f.v k (dirs.take k) f.done false).2.2
          (refOuts g f.v k (dirs.take k) f.done false).1
          (fun h => (refOuts_hv_false _ _ h).2)
        refine ⟨none, A'' ++ A', B' ++ B'', fun ps => ?_⟩
        obtain ⟨h1, h2⟩ := ih ps
        rcases hx : openBlockP g dirs (k + 1) fs ps with ⟨r₀, ps₀⟩
        rw [hx] at h1 h2
        simp only at h1 h2
        subst h1
        simp only [openBlockP, hx, hk, ↓reduceIte]
        refine ⟨by first | trivial | rfl, ?_⟩
        rw [ho, hf ps₀, h2]
        simp only [List.append_assoc]

theorem openBlock_ctx (g : Graph) (fs : List PathFrame) (dirs : List Bool) (hT : TreeFrames fs) :
    ∃ r A B, ∀ ps, openBlock g fs dirs ps = ⟨r, A ++ stNest ps ++ B⟩ := by
  obtain ⟨r, A, B, h⟩ := openBlockP_ctx g dirs fs 0 hT
  exact ⟨r, A, B, fun ps => by rw [openBlock, (h ps).1, (h ps).2]⟩

theorem stNest_snoc (ps : List StPiece) (p : StPiece) :
    stNest (ps ++ [p]) = stNestL [p] ++ stNest ps ++ stNestR [p] := by
  simp [stNest, stNestL_append, stNestR_append, stNestL, stNestR]

end Spqr

namespace Spqr

/-- Well-formedness of a frame of the open path relative to the vertices `verts` / edges `edges`
still to be walked below it: a tree edge, distinct and in-range finished vertices, disjoint
from what is still to come. -/
structure FrameOK (g : Graph) (f : PathFrame) (verts edges : List Nat) : Prop where
  tree : ∃ e cls child, f.o = .tree e cls child
  vnodup : (f.v :: DfsOut.vertsList f.done).Nodup
  enodup : (DfsOut.edgesList f.done).Nodup
  vlt : ∀ w ∈ f.v :: DfsOut.vertsList f.done, w < g.nv
  vdisj : ∀ w ∈ f.v :: DfsOut.vertsList f.done, w ∉ verts
  edisj : ∀ e ∈ f.o.e :: DfsOut.edgesList f.done, e ∉ edges

theorem FrameOK.mono {g : Graph} {f : PathFrame} {verts edges verts' edges' : List Nat}
    (h : FrameOK g f verts edges) (hv : ∀ w ∈ verts', w ∈ verts) (he : ∀ e ∈ edges', e ∈ edges) :
    FrameOK g f verts' edges' :=
  ⟨h.tree, h.vnodup, h.enodup, h.vlt, fun w hw hw' => h.vdisj w hw (hv w hw'),
    fun e he' he'' => h.edisj e he' (he e he'')⟩

theorem TreeFrames.of_ok {g : Graph} {fs : List PathFrame} {verts edges : List Nat}
    (h : ∀ f ∈ fs, FrameOK g f verts edges) : TreeFrames fs := fun f hf => (h f hf).tree

end Spqr
