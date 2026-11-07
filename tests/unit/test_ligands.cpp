#include "mini_test.h"
#include "leak_probe.h"
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

/* ---- F5: whitespace/tab separated ATOM records ---- */
static void checkSameAsSelf(const char* fixture) {
    Ligands fixed("fixtures/self.mol2", false);
    Ligands loose(fixture, false);
    std::vector<Ligand> a = fixed.getLigands();
    std::vector<Ligand> b = loose.getLigands();
    CHECK(loose.size() == 1);
    CHECK(b[0].size() == a[0].size());
    for (int i = 0; i < (int)a[0].size() && i < (int)b[0].size(); i++) {
        CHECK_CLOSE(b[0].getIthPoints(i)[0], a[0].getIthPoints(i)[0], 1e-6);
        CHECK_CLOSE(b[0].getIthPoints(i)[1], a[0].getIthPoints(i)[1], 1e-6);
        CHECK_CLOSE(b[0].getIthPoints(i)[2], a[0].getIthPoints(i)[2], 1e-6);
        CHECK(b[0].getAtomType(i) == a[0].getAtomType(i));
        CHECK(b[0].getAtomName(i) == a[0].getAtomName(i));
    }
}

TEST(ligands_space_delimited_matches_fixed_columns) {
    checkSameAsSelf("fixtures/space_delimited.mol2");
}

TEST(ligands_tab_delimited_matches_fixed_columns) {
    checkSameAsSelf("fixtures/tab_delimited.mol2");
}

TEST(ligands_unknown_extension_falls_back_to_mol2) {
    /* fixtures/self.pdb is a mol2 file with a misleading extension; the
     * factory warns on stderr and still parses it as mol2. */
    Ligands ligs("fixtures/self.pdb", false);
    CHECK(ligs.size() == 1);
    CHECK(ligs[0].size() == 6);
    CHECK(ligs[0].getLigname() == "SELF");
}

/* ---- B29: original atom names are preserved ---- */
TEST(ligands_atom_names_preserved) {
    Ligands ligs("fixtures/atom_name_align.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    CHECK(all.size() == 1);
    CHECK(all[0].size() == 5);
    CHECK(all[0].getAtomName(0) == "Cl1");
    CHECK(all[0].getAtomName(1) == "1H");
    CHECK(all[0].getAtomName(2) == "C10");
    CHECK(all[0].getAtomName(3) == "ABCDE");
    CHECK(all[0].getAtomName(4) == "");
    CHECK(all[0].getAtomName(-1) == "");
    CHECK(all[0].getAtomName(99) == "");
}

TEST(ligands_atom_names_survive_copy_and_assign) {
    Ligands ligs("fixtures/self.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand copy(all[0]);
    CHECK(copy.getAtomName(0) == "C1");
    CHECK(copy.getAtomName(5) == "C6");

    Ligand assigned;
    assigned = all[0];
    CHECK(assigned.getAtomName(2) == "C3");
    CHECK(assigned.getAtomName(2) == all[0].getAtomName(2));
}

TEST(ligands_atom_names_dropped_with_hydrogens) {
    Ligands withH("../example/xtal-lig.mol2", true);
    Ligands noH("../example/xtal-lig.mol2", false);
    std::vector<Ligand> wh = withH.getLigands();
    std::vector<Ligand> nh = noH.getLigands();
    CHECK(wh[0].size() == 59);
    CHECK(nh[0].size() == 29);
    /* the kept names must follow the kept atoms (no index drift) */
    int kept = 0;
    for (int i = 0; i < (int)wh[0].size(); i++) {
        if (wh[0].getAtomType(i) != "H") {
            CHECK(nh[0].getAtomName(kept) == wh[0].getAtomName(i));
            kept++;
        }
    }
    CHECK(kept == 29);
    CHECK(nh[0].getAtomName(28) == "C21");
}

/* ---- B29: PDB name alignment and element column ---- */
static std::vector<std::string> pdbOf(const char* fixture) {
    Ligands ligs(fixture, false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    return L.format2PDB("QUE");
}

TEST(pdb_atom_name_alignment_rules) {
    std::vector<std::string> rows = pdbOf("fixtures/atom_name_align.mol2");
    CHECK(rows.size() == 6);                       /* 5 atoms + TER */
    CHECK(rows[0].substr(12, 4) == "Cl1 ");        /* two-letter element */
    CHECK(rows[1].substr(12, 4) == "1H  ");        /* leading digit */
    CHECK(rows[2].substr(12, 4) == " C10");        /* single letter, right */
    CHECK(rows[3].substr(12, 4) == "ABCD");        /* truncated to 4 chars */
    CHECK(rows[4].substr(12, 4) == " O  ");        /* empty -> element */
    for (int i = 0; i < 5; i++) {
        CHECK(rows[i].substr(17, 3) == "QUE");     /* residue not displaced */
        CHECK(rows[i].size() == 78);
    }
    CHECK(rows[0].substr(76, 2) == "CL");          /* Cl -> 77-78 */
    CHECK(rows[1].substr(76, 2) == " N");          /* N  -> 77-78, right */
    CHECK(rows[2].substr(76, 2) == " C");
    CHECK(rows[3].substr(76, 2) == " C");
    CHECK(rows[4].substr(76, 2) == " O");
}

TEST(pdb_all_thin_overloads_agree_on_names) {
    Ligands ligs("fixtures/self.mol2", false);
    std::vector<Ligand> all = ligs.getLigands();
    Ligand L = all[0];
    std::vector<std::string> a = L.format2PDB();
    std::vector<std::string> b = L.format2PDB("QUE");
    std::vector<std::string> c = L.format2PDB("QUE", 1, 88888);
    CHECK(a.size() == b.size() && b.size() == c.size());
    for (int i = 0; i + 1 < (int)c.size(); i++) { /* last line is TER */
        CHECK(b[i] == c[i]);                       /* full line */
        /* the no-arg overload differs only in residue/serial columns */
        CHECK(a[i].substr(12, 4) == c[i].substr(12, 4));
        CHECK(a[i].substr(76, 2) == c[i].substr(76, 2));
    }
    CHECK(a.back() == "TER");
    /* no-arg overload keeps its historical JUN residue */
    CHECK(a[0].substr(17, 3) == "JUN");
    CHECK(b[0].substr(17, 3) == "QUE");
}

/* ---- F6: a zero-atom ligand must not rebuild its lazy cache every call ---- */
TEST(f6_empty_ligand_builds_once) {
    Ligands ligs("fixtures/empty_ligand.mol2", false);
    std::vector<Ligand> objs;
    objs.reserve(100);
    for (int k = 0; k < 100; ++k) objs.push_back(ligs[0]);
    { std::vector<std::vector<BondType> > warm = objs[0].getAtomBondTypes(); (void)warm; }
    LeakDelta d;
    for (int k = 0; k < 100; ++k) {
        const std::vector<std::vector<BondType> >& abt = objs[k].getAtomBondTypes();
        CHECK(abt.size() == 0);
    }
    CHECK(d.delta() == 0);
}

/* ---- F6: the getters hand out the stored vectors, not copies ---- */
TEST(f6_const_ref_no_copy) {
    Ligands ligs("fixtures/isolated_atom.mol2", false);
    const Ligand& L = ligs[0];
    /* read-only getters are const-qualified and return references */
    CHECK(&L.getPoints() == &L.getPoints());
    CHECK(&L.getAtomTypes() == &L.getAtomTypes());
    CHECK(&L.getBondInfoes() == &L.getBondInfoes());
    /* getAtomBondTypes() builds its cache lazily, so it needs a non-const
     * object (unchanged by this ticket); the reference it returns is stable */
    Ligand M = L;
    const std::vector<std::vector<BondType> >& a = M.getAtomBondTypes();
    const std::vector<std::vector<BondType> >& b = M.getAtomBondTypes();
    CHECK(&a == &b);
    /* Ligands::operator[] returns the stored element, not a fresh copy: two
     * such addresses compared in one full-expression would differ otherwise. */
    CHECK(&ligs[0] == &ligs[0]);
}

/* ---- F6: atomBondTypes.size() == points.size() even with isolated atoms ---- */
TEST(f6_size_invariant) {
    const char* fixtures[3] = {"fixtures/isolated_atom.mol2",
                               "fixtures/empty_ligand.mol2",
                               "fixtures/self.mol2"};
    for (int f = 0; f < 3; ++f) {
        Ligands ligs(fixtures[f], false);
        Ligand L = ligs[0];
        const std::vector<std::vector<BondType> >& abt = L.getAtomBondTypes();
        CHECK((int)abt.size() == L.size());
        CHECK(abt.size() == L.getPoints().size());
    }

    /* isolated_atom.mol2: atoms 1-3 form a chain, atom 4 has no bond */
    {
        Ligands iso("fixtures/isolated_atom.mol2", false);
        Ligand L = iso[0];
        const std::vector<std::vector<BondType> >& abt = L.getAtomBondTypes();
        CHECK(abt.size() == 4);
        CHECK(abt[0].size() == 1);
        CHECK(abt[1].size() == 2);
        CHECK(abt[3].empty());
    }
}
