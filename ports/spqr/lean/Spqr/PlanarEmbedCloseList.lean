import Spqr.PlanarEmbedClose

namespace Spqr.PlanarSpqrTree

def closeList (js : List Nat) (s : EmbedState) : EmbedState :=
  js.foldl (fun s j => ((closeOuter j).run s).2) s

theorem closeList_nil (s : EmbedState) : closeList [] s = s := rfl

theorem closeList_cons (j : Nat) (js : List Nat) (s : EmbedState) :
    closeList (j :: js) s = closeList js ((closeOuter j).run s).2 := rfl

theorem closeList_outerE (js : List Nat) (s : EmbedState) : (closeList js s).outerE = s.outerE := by
  induction js generalizing s with
  | nil => rfl
  | cons j js ih => rw [closeList_cons, ih, closeOuter_outerE]

theorem closeList_rotAdj_size (js : List Nat) (s : EmbedState) :
    (closeList js s).rotAdj.size = s.rotAdj.size := by
  induction js generalizing s with
  | nil => rfl
  | cons j js ih => rw [closeList_cons, ih, closeOuter_rotAdj_size]

theorem closeOuter_get_congr (j q : Nat) (s s' : EmbedState)
    (ho : s'.outerE = s.outerE) (hs : s'.rotAdj.size = s.rotAdj.size)
    (hq : s'.rotAdj[q]? = s.rotAdj[q]?) :
    ((closeOuter j).run s').2.rotAdj[q]? = ((closeOuter j).run s).2.rotAdj[q]? := by
  rw [closeOuter_run, closeOuter_run, ho]
  cases s.outerE[j]![0]! with
  | none => simpa only [link_none_left] using hq
  | some a =>
    cases s.outerE[j]![1]! with
    | none => simpa only [link_none_right] using hq
    | some b =>
      rw [link_some_run, link_some_run]
      simp only [Array.set!, Array.getElem?_setIfInBounds, Array.size_setIfInBounds, hs, hq]

theorem closeList_frame (js : List Nat) (s : EmbedState) (q : Nat)
    (h : ∀ j ∈ js, ¬s.exposedAt j q) : (closeList js s).rotAdj[q]? = s.rotAdj[q]? := by
  induction js generalizing s with
  | nil => rfl
  | cons j js ih =>
    rw [closeList_cons, ih, closeOuter_frame j s q (h j (by simp))]
    intro k hk
    simpa only [EmbedState.exposedAt, closeOuter_outerE] using h k (by simp [hk])

theorem closeList_get (js : List Nat) (hnd : js.Nodup) (j : Nat) (hj : j ∈ js)
    (s : EmbedState) (q : Nat)
    (hothers : ∀ k ∈ js, k ≠ j → ¬s.exposedAt k q) :
    (closeList js s).rotAdj[q]? = ((closeOuter j).run s).2.rotAdj[q]? := by
  induction js generalizing s with
  | nil => simp at hj
  | cons k js ih =>
    rw [List.nodup_cons] at hnd
    rw [closeList_cons]
    rcases List.mem_cons.1 hj with hjk | hjt
    · subst k
      apply closeList_frame
      intro k hk
      have hkj : k ≠ j := by intro heq; subst k; exact hnd.1 hk
      simpa only [EmbedState.exposedAt, closeOuter_outerE] using hothers k (by simp [hk]) hkj
    · rw [ih hnd.2 hjt]
      · apply closeOuter_get_congr
        · exact closeOuter_outerE k s
        · exact closeOuter_rotAdj_size k s
        · exact closeOuter_frame k s q (hothers k (by simp) (by intro heq; subst k; exact hnd.1 hjt))
      · intro l hl hlj
        simpa only [EmbedState.exposedAt, closeOuter_outerE] using hothers l (by simp [hl]) hlj

theorem forM_closeOuter_run (js : List Nat) (s : EmbedState) :
    ((js.forM closeOuter).run s).2 = closeList js s := by
  induction js generalizing s with
  | nil => rfl
  | cons j js ih =>
    change ((js.forM closeOuter).run ((closeOuter j).run s).2).2 = _
    exact ih _

theorem embedItem_F (t : PlanarSpqrTree) (i : Nat) (hty : t.types[i]! = .F) (s : EmbedState) :
    ((t.embedItem i).run s).2 = closeList (t.children i) s := by
  have idbind {α β : Type} (a : Id α) (f : α → Id β) : (a >>= f) = f a := rfl
  have idpure {α : Type} (a : α) : (pure a : Id α) = a := rfl
  have idmap {α β : Type} (a : Id α) (f : α → β) : (f <$> a) = f a := rfl
  unfold embedItem
  rw [hty]
  simp only
  generalize t.children i = js
  induction js generalizing s with
  | nil => rfl
  | cons j js ih =>
    cases ha : s.outerE[j]![0]! <;> cases hb : s.outerE[j]![1]! <;>
      simpa [List.forIn_cons, closeList_cons, closeOuter_run, outer_run,
        ha, hb, StateT.run_bind, link_none_left, link_none_right, link_some_run,
        idbind, idpure, idmap] using
        ih ((closeOuter j).run s).2

end Spqr.PlanarSpqrTree
