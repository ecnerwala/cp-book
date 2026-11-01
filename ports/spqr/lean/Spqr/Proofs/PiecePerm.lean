import Spqr.Proofs.PieceInsert
import Spqr.Proofs.PlanarReindex

/-!
# Reordering the edge list of a `Capped` piece

`Piece.Capped` is stated for a particular order of `ves` (through `Piece.loc`). The node fold builds
its certificate in the order the children are glued, while `pieceBelow` lists the children in tree
order; `Capped.perm` transports the certificate along a permutation of `ves` (both lists nodup) by
reindexing the planar witness.
-/

namespace Spqr

namespace Piece

variable {P : Piece} {q l : Nat}

/-- The global quarter-edge of local quarter-edge `l`. -/
def glob (P : Piece) (l : Nat) : Nat := 4 * (P.ves[l / 4]?.getD 0) + l % 4

theorem glob_loc (h : P.loc q = some l) : P.glob l = q := by
  obtain ⟨k, hk, hkq, rfl⟩ := loc_data h
  have e1 : (4 * k + q % 4) / 4 = k := by omega
  have e2 : (4 * k + q % 4) % 4 = q % 4 := by omega
  simp only [glob, e1, e2, hkq, Option.getD_some]
  unfold QE.edge; omega

theorem loc_glob (hnd : P.ves.Nodup) (hl : l < 4 * P.ves.length) : P.loc (P.glob l) = some l := by
  have hk : l / 4 < P.ves.length := by omega
  obtain ⟨x, hx⟩ : ∃ x, P.ves[l / 4]? = some x := ⟨_, List.getElem?_eq_getElem hk⟩
  have hx' : P.ves[l / 4] = x := by simpa [List.getElem?_eq_getElem hk] using hx
  have hg : P.glob l = 4 * x + l % 4 := by simp [glob, hx]
  have he : QE.edge (P.glob l) = x := by unfold QE.edge; omega
  have hfind : P.ves.findIdx? (· == QE.edge (P.glob l)) = some (l / 4) := by
    rw [List.findIdx?_eq_some_iff_getElem]
    refine ⟨hk, by rw [he, hx']; simp, ?_⟩
    intro j hj
    rw [he]
    simp only [beq_iff_eq]
    intro heq
    have := (List.Nodup.getElem_inj_iff hnd).1 (heq.trans hx'.symm)
    omega
  show (P.ves.findIdx? (· == QE.edge (P.glob l))).map (fun k => 4 * k + P.glob l % 4) = some l
  rw [hfind, hg]
  simp only [Option.map_some, Option.some.injEq]
  omega

theorem glob_xor {c : Nat} (hc : c < 4) : P.glob (l ^^^ c) = P.glob l ^^^ c := by
  simp only [glob, xor_div4 l hc, xor_mod4 l hc]
  exact (mul4_add_xor _ (l % 4) (Nat.mod_lt _ (by decide)) hc).symm

theorem glob_mod_two : P.glob l % 2 = l % 2 := by
  unfold glob; omega

end Piece

/-- `faceStep` is natural under a quarter-edge bijection compatible with `get` and `^^^ 3`. -/
theorem sameFaceOrbit_of_map {ρ σ : RotationSystem} (φ : Nat → Nat)
    (hget : ∀ q, q < σ.size → ρ.get (φ q) = (σ.get q).map φ)
    (hbound : ∀ q, q < σ.size → ∀ r ∈ σ.get q, r < σ.size)
    (hinj : ∀ q r, q < σ.size → r < σ.size → φ q = φ r → q = r)
    (hxor : ∀ q, q < σ.size → φ (q ^^^ 3) = φ q ^^^ 3)
    (hs4 : ∃ m, σ.size = 4 * m)
    {a b : Nat} (ha : a < σ.size) (hb : b < σ.size) (h : ρ.SameFaceOrbit (φ a) (φ b)) :
    σ.SameFaceOrbit a b := by
  obtain ⟨m, hm⟩ := hs4
  have key : ∀ x y, Relation.ReflTransGen (fun a b => ρ.faceStep a = some b) x y →
      ∀ a, a < σ.size → φ a = x → ∃ b', b' < σ.size ∧ φ b' = y ∧ σ.SameFaceOrbit a b' := by
    intro x y hxy
    induction hxy with
    | refl => intro a ha hx; exact ⟨a, ha, hx, Relation.ReflTransGen.refl⟩
    | tail _ hst ih =>
      intro a ha hx
      obtain ⟨b', hb', rfl, hab⟩ := ih a ha hx
      have hb3 : b' ^^^ 3 < σ.size := by rw [hm] at hb' ⊢; exact xor_lt_mul4 hb' (by decide)
      have hstep := hst
      unfold RotationSystem.faceStep at hstep
      rw [show QE.across (φ b') = φ b' ^^^ 3 from rfl, ← hxor b' hb', hget _ hb3] at hstep
      obtain ⟨c, hc, rfl⟩ := Option.map_eq_some_iff.1 hstep
      refine ⟨c, hbound _ hb3 c hc, rfl, Relation.ReflTransGen.tail hab ?_⟩
      unfold RotationSystem.faceStep
      exact hc
  obtain ⟨b', hb', hφ, hab⟩ := key _ _ h a ha rfl
  rwa [hinj b' b hb' hb hφ] at hab

namespace Piece

variable {P : Piece} {E : List Nat}

theorem loc_of_loc_perm (hperm : P.ves.Perm E) (hnd : P.ves.Nodup) {q l : Nat}
    (h : ({P with ves := E} : Piece).loc q = some l) :
    P.loc q = some ((P.loc (({P with ves := E} : Piece).glob l)).getD 0) := by
  rw [glob_loc h]
  have hm : P.Mem q := by
    have := mem_of_loc h
    unfold Mem at this ⊢
    exact hperm.mem_iff.2 this
  obtain ⟨l', hl'⟩ := loc_exists hm
  rw [hl']; rfl

/-- `Capped` is invariant under reordering `ves`. -/
theorem Capped.perm {A : Array (Option Nat)} {ρ : RotationSystem} {c0 c1 c2 c3 u v : Nat}
    (h : P.Capped A ρ c0 c1 c2 c3 u v) (hperm : P.ves.Perm E) (hnd : P.ves.Nodup) :
    ∃ σ, ({P with ves := E} : Piece).Capped A σ c0 c1 c2 c3 u v := by
  set Q : Piece := {P with ves := E} with hQ
  have hndE : E.Nodup := hperm.nodup_iff.1 hnd
  have hlen : P.ves.length = E.length := hperm.length_eq
  have hρs : ρ.size = 4 * P.ves.length := by rw [h.planar.size]; simp [es]
  let φ : Nat → Nat := fun l => (P.loc (Q.glob l)).getD 0
  let ψ : Nat → Nat := fun r => (Q.loc (P.glob r)).getD 0
  have hmemQP : ∀ q, Q.Mem q ↔ P.Mem q := fun q => by
    unfold Mem; exact hperm.mem_iff.symm
  -- Q.loc q = some l  →  P.loc q = some (φ l)
  have hQP : ∀ {q l}, Q.loc q = some l → P.loc q = some (φ l) := fun hl =>
    loc_of_loc_perm hperm hnd hl
  -- P.loc q = some r  →  Q.loc q = some (ψ r)
  have hPQ : ∀ {q r}, P.loc q = some r → Q.loc q = some (ψ r) := by
    intro q r hr
    show Q.loc q = some ((Q.loc (P.glob r)).getD 0)
    rw [glob_loc hr]
    obtain ⟨l', hl'⟩ := loc_exists ((hmemQP q).2 (mem_of_loc hr))
    rw [hl']; rfl
  have hφlt : ∀ l, l < 4 * E.length → φ l < ρ.size := by
    intro l hl
    have := hQP (loc_glob (P := Q) hndE hl)
    rw [hρs]; exact loc_lt this
  have hψlt : ∀ r, r < ρ.size → ψ r < 4 * E.length := by
    intro r hr
    rw [hρs] at hr
    have := hPQ (loc_glob hnd hr)
    exact loc_lt this
  have hφψ : ∀ r, r < ρ.size → φ (ψ r) = r := by
    intro r hr
    rw [hρs] at hr
    have h1 := hPQ (loc_glob hnd hr)
    have h2 := hQP h1
    rw [loc_glob hnd hr] at h2
    exact (Option.some.inj h2).symm
  have hψφ : ∀ l, l < 4 * E.length → ψ (φ l) = l := by
    intro l hl
    have h1 := hQP (loc_glob (P := Q) hndE hl)
    have h2 := hPQ h1
    rw [loc_glob (P := Q) hndE hl] at h2
    exact (Option.some.inj h2).symm
  let σ : RotationSystem := ⟨Array.ofFn (n := 4 * E.length) fun q => (ρ.get (φ q)).map ψ⟩
  have hσs : σ.size = 4 * E.length := by simp [σ, RotationSystem.size]
  have hσget : ∀ q, q < σ.size → σ.get q = (ρ.get (φ q)).map ψ := by
    intro q hq
    rw [hσs] at hq
    simp [σ, RotationSystem.get, Array.getElem?_ofFn, hq]
  have hbound : ∀ q, q < σ.size → ∀ r ∈ σ.get q, r < σ.size := by
    intro q hq r hr
    rw [hσget q hq] at hr
    obtain ⟨r', hr', rfl⟩ := Option.map_eq_some_iff.1 hr
    have := (h.planar.involution _ (hφlt q (hσs ▸ hq)) r' hr').1
    rw [hσs]; exact hψlt r' this
  have hget : ∀ q, q < σ.size → ρ.get (φ q) = (σ.get q).map φ := by
    intro q hq
    rw [hσget q hq, Option.map_map]
    cases hr : ρ.get (φ q) with
    | none => rfl
    | some r =>
      have := (h.planar.involution _ (hφlt q (hσs ▸ hq)) r hr).1
      simp [Function.comp, hφψ r this]
  have hinj : ∀ q r, q < σ.size → r < σ.size → φ q = φ r → q = r := by
    intro q r hq hr he
    have := congrArg ψ he
    rwa [hψφ q (hσs ▸ hq), hψφ r (hσs ▸ hr)] at this
  have hxor : ∀ q, q < σ.size → ∀ c, c < 4 → φ (q ^^^ c) = φ q ^^^ c := by
    intro q hq c hc
    have h1 := hQP (loc_glob (P := Q) hndE (hσs ▸ hq))
    have h2 := loc_xor h1 hc
    show (P.loc (Q.glob (q ^^^ c))).getD 0 = φ q ^^^ c
    rw [glob_xor hc, h2]; rfl
  have hdir : ∀ q, q < σ.size → QE.dir (φ q) = QE.dir q := by
    intro q hq
    have h1 := hQP (loc_glob (P := Q) hndE (hσs ▸ hq))
    have := loc_mod_two h1
    unfold QE.dir
    rw [this, glob_mod_two]
  have hvert : ∀ q, q < σ.size → QE.vert P.es (φ q) = QE.vert Q.es q := by
    intro q hq
    have h0 := loc_glob (P := Q) hndE (hσs ▸ hq)
    have h1 := hQP h0
    rw [vert_loc h1, vert_loc h0]
  have hmem : ∀ p, p ∈ P.es ↔ p ∈ Q.es := fun p => (hperm.map P.ends).mem_iff
  have hplanar : IsPlanarEmbedding Q.es Q.nVerts σ :=
    h.planar.reindex φ (by simp [hσs, Q, es]) (by simp [Q, es, hlen]) hmem (fun q hq => hφlt q (hσs ▸ hq)) hinj hget
      hbound hdir hvert hxor
  have hσt : σ.Total := hplanar.total
  have hσi : σ.Involution := hplanar.involution
  -- transport of `get`
  have hgetT : ∀ {lq lr}, ρ.get lq = some lr → lq < ρ.size → σ.get (ψ lq) = some (ψ lr) := by
    intro lq lr hlr hlq
    rw [hσget _ (hσs ▸ hψlt lq hlq), hφψ lq hlq, hlr]; rfl
  have hvertT : ∀ {l}, l < ρ.size → QE.vert Q.es (ψ l) = QE.vert P.es l := by
    intro l hl
    rw [← hvert _ (hσs ▸ hψlt l hl), hφψ l hl]
  refine ⟨σ, hplanar, ?_, ?_, ?_, ?_, ?_, h.dir0, h.dir1, h.dir2, h.dir3⟩
  · intro q r hq hr
    obtain ⟨lq, lr, hlq, hlr, hg⟩ := h.agrees q r ((hmemQP q).1 hq) hr
    exact ⟨ψ lq, ψ lr, hPQ hlq, hPQ hlr, hgetT hg (hρs ▸ loc_lt hlq)⟩
  · intro q hq
    exact h.unset q ((hmemQP q).1 hq)
  · obtain ⟨l0, l1, h0, h1, hg, hv⟩ := h.pair0
    exact ⟨ψ l0, ψ l1, hPQ h0, hPQ h1, hgetT hg (hρs ▸ loc_lt h0),
      by rw [hvertT (hρs ▸ loc_lt h0)]; exact hv⟩
  · obtain ⟨l2, l3, h2, h3, hg, hv⟩ := h.pair2
    exact ⟨ψ l2, ψ l3, hPQ h2, hPQ h3, hgetT hg (hρs ▸ loc_lt h2),
      by rw [hvertT (hρs ▸ loc_lt h2)]; exact hv⟩
  · obtain ⟨l0, l2, h0, h2, hf⟩ := h.face
    have hl0 := hρs ▸ loc_lt h0
    have hl2 := hρs ▸ loc_lt h2
    refine ⟨ψ l0, ψ l2, hPQ h0, hPQ h2, ?_⟩
    refine sameFaceOrbit_of_map φ hget hbound hinj (fun q hq => hxor q hq 3 (by decide))
      ⟨E.length, hσs⟩ (hσs ▸ hψlt l0 hl0) (hσs ▸ hψlt l2 hl2) ?_
    rw [hφψ l0 hl0, hφψ l2 hl2]; exact hf

end Piece

end Spqr
