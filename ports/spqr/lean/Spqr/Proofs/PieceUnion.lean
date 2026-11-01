import Spqr.Proofs.PieceAppend
import Spqr.PlanarEmbedLeaf

namespace Spqr.Piece

theorem union_embeddings (P : Piece) (js : List Nat) (edges : Nat → List Nat)
    (a : Array (Option Nat))
    (hpieces : ∀ j ∈ js, ∃ ρ, IsPlanarEmbedding ({P with ves := edges j} : Piece).es P.nVerts ρ ∧
      ({P with ves := edges j} : Piece).Agrees a ρ)
    (hdis : js.Pairwise fun j k => List.Disjoint (edges j) (edges k) ∧
      ∀ v, HasEdge ({P with ves := edges j} : Piece).es v →
        HasEdge ({P with ves := edges k} : Piece).es v → False) :
    ∃ ρ, IsPlanarEmbedding ({P with ves := js.flatMap edges} : Piece).es P.nVerts ρ ∧
      ({P with ves := js.flatMap edges} : Piece).Agrees a ρ := by
  induction js with
  | nil =>
    refine ⟨⟨#[]⟩, isPlanarEmbedding_nil _, ?_⟩
    intro q r hq
    simp [Mem] at hq
  | cons j js ih =>
    obtain ⟨ρ₁, hρ₁, ha₁⟩ := hpieces j (by simp)
    obtain ⟨hp, ht⟩ := List.pairwise_cons.1 hdis
    obtain ⟨ρ₂, hρ₂, ha₂⟩ := ih (fun k hk => hpieces k (by simp [hk])) ht
    have hed : List.Disjoint (edges j) (js.flatMap edges) := by
      rw [List.disjoint_left]
      intro e he hf
      obtain ⟨k, hk, he'⟩ := List.mem_flatMap.1 hf
      exact List.disjoint_left.1 (hp k hk).1 he he'
    have hvd : ∀ v, HasEdge ({P with ves := edges j} : Piece).es v →
        HasEdge ({P with ves := js.flatMap edges} : Piece).es v → False := by
      rintro v hv ⟨p, hp', hv'⟩
      obtain ⟨e, he, rfl⟩ := List.mem_map.1 hp'
      obtain ⟨k, hk, he'⟩ := List.mem_flatMap.1 he
      exact (hp k hk).2 v hv ⟨P.ends e, List.mem_map.2 ⟨e, he', rfl⟩, hv'⟩
    refine ⟨ρ₁.union ρ₂, ?_, ?_⟩
    · simpa only [List.flatMap_cons, es, List.map_append] using hρ₁.append hρ₂ hvd
    · exact agrees_append hed (by simpa [es] using hρ₁.size) ha₁ ha₂

end Spqr.Piece
