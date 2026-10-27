#include "mini_test.h"
#include "Ligands.h"
#include <string>

/* Fixtures are parsed exactly as the program parses real mol2 files
 * (fixed columns, atom type stripped to element symbols, H handling). */

TEST(ligands_loads_single_fixture) {
    Ligands ligs("fixtures/self.mol2", false);
    CHECK(ligs.size() == 1);
    CHECK(ligs[0].size() == 6);
    CHECK(ligs[0].getLigname() == "SELF");
}

TEST(ligands_atom_type_reduced_to_element) {
    Ligands ligs("fixtures/weird_cols.mol2", false);
    CHECK(ligs.size() == 1);
    /* getLigands() returns by value -- hold the copy, do not bind a reference */
    std::vector<Ligand> all = ligs.getLigands();
    Ligand& L = all[0];
    /* C.ar -> C, N.am -> N, O.co2 -> O, Cl -> CL, Br -> BR (uppercased) */
    CHECK(L.getAtomType(0) == "C");
    CHECK(L.getAtomType(1) == "N");
    CHECK(L.getAtomType(2) == "O");
    CHECK(L.getAtomType(3) == "CL");
    CHECK(L.getAtomType(4) == "BR");
}

TEST(ligands_bond_count) {
    Ligands ligs("fixtures/self.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    CHECK(all[0].getBondInfoes().size() == 6);
}

TEST(ligands_multiligand_file) {
    /* example/ace_ligands.mol2 has 49 ligands */
    Ligands ligs("../example/ace_ligands.mol2", false);
    CHECK(ligs.size() == 49);
    /* first ligand is ZINC01535869 with 29 heavy atoms (H dropped) */
    CHECK(ligs[0].getLigname() == "ZINC01535869");
}

TEST(ligands_hydrogens_dropped_by_default) {
    /* xtal-lig.mol2: 59 atoms incl. H -> 29 heavy atoms without H */
    Ligands noH("../example/xtal-lig.mol2", false);
    CHECK(noH[0].size() == 29);
    Ligands withH("../example/xtal-lig.mol2", true);
    CHECK(withH[0].size() == 59);
}

/* ---- B24: an empty ATOM section parses to a zero-atom ligand ---- */
TEST(empty_ligand_from_parse_has_zero_atoms) {
    /* Parse only -- no alignment.  Before the fix this already did not crash,
     * so it is a plain TEST regression lock, not an XFAIL. */
    Ligands ligs("fixtures/empty_ligand.mol2", false);
    CHECK(ligs.size() == 1);
    CHECK(ligs[0].size() == 0);
    CHECK(ligs[0].getLigname() == "EMPTYLIG");
}
