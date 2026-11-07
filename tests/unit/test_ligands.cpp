#include "mini_test.h"
#include "leak_probe.h"
#include "Ligands.h"
#include <string>
#include <vector>

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
    /* getLigands() returns const& -- the copy below is a deliberate snapshot */
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
    /* Ligands::operator[] must return the stored element, not a fresh copy:
     * binding the same element twice must yield one address.  (A by-value
     * operator[] binds two distinct temporaries and turns this red.) */
    const Ligand& first = ligs[0];
    const Ligand& second = ligs[0];
    CHECK(&first == &second);
    /* and a getter reached through operator[] is not a reference into a
     * temporary: it must be the vector owned by the stored element */
    CHECK(&ligs[0].getPoints() == &first.getPoints());
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

/* ---- F11: setPoints must drop the caches it invalidated ---- */
TEST(f11_setpoints_invalidates_caches) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligand L = ligs[0];

    /* warm both lazy caches (bondDisMtx, btMtx and atomBondTypes) */
    CHECK(L.getBondDistMtx() != NULL);
    CHECK((int)L.getAtomBondTypes().size() == L.size());
    long warmed = liveAllocs();

    /* a rigid shift: same atom count, every bond length is preserved */
    std::vector<Point3D> shifted;
    for (int i = 0; i < (int)L.size(); ++i) {
        Point3D p = L.getIthPoints(i);
        shifted.push_back(Point3D(p[0] + 10.0, p[1] - 3.0, p[2] + 7.0));
    }
    L.setPoints(shifted);

    /* the matrices were released (7 allocations each) and are NULL again,
     * so the getters rebuild them lazily from the new geometry.  Without
     * resetCaches() the stale matrices stay allocated and this is red. */
    long after = liveAllocs();
    CHECK(after + 5 <= warmed);

    double** reborn = L.getBondDistMtx();
    CHECK(reborn != NULL);
    CHECK(reborn[0][0] == -1.0);               /* rebuilt from scratch */
    CHECK(reborn[0][1] > 0.0);                 /* ring bond 0-1 survived */

    /* atomBondTypes was cleared too and is rebuilt on demand */
    CHECK((int)L.getAtomBondTypes().size() == L.size());
}

/* ---- F15: bondDis follows the coordinates on a non-rigid setPoints ---- */
TEST(f15_recompute_nonrigid) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligand L = ligs[0];

    std::vector<double> before;
    const std::vector<BondInfo>& bi = L.getBondInfoes();
    CHECK(bi.size() == 6);
    for (int i = 0; i < (int)bi.size(); ++i) before.push_back(bi[i].bondDis);
    CHECK(before[0] > 0.0);

    /* scaling doubles every bond length: the one transform that a rigid
     * check cannot catch.  (Before the fix setPoints left the old bondDis,
     * so this is red: the values stay at the unscaled lengths.) */
    std::vector<Point3D> scaled;
    for (int i = 0; i < (int)L.size(); ++i) {
        Point3D p = L.getIthPoints(i);
        scaled.push_back(Point3D(p[0] * 2.0, p[1] * 2.0, p[2] * 2.0));
    }
    L.setPoints(scaled);

    const std::vector<BondInfo>& after = L.getBondInfoes();
    CHECK(after.size() == bi.size());
    for (int i = 0; i < (int)after.size(); ++i) {
        CHECK_CLOSE(after[i].bondDis, before[i] * 2.0, 1e-9);
        CHECK_CLOSE(after[i].bondDisPow2, before[i] * before[i] * 4.0, 1e-9);
    }
    /* the distance matrix reads the refreshed bondDis */
    double** mtx = L.getBondDistMtx();
    CHECK_CLOSE(mtx[0][1], before[0] * 2.0, 1e-9);
}

/* ---- F15: a rigid transform leaves every bond length alone ---- */
TEST(f15_rigid_unchanged) {
    Ligands ligs("fixtures/self.mol2", false);
    Ligand L = ligs[0];

    std::vector<double> before;
    const std::vector<BondInfo>& bi = L.getBondInfoes();
    for (int i = 0; i < (int)bi.size(); ++i) before.push_back(bi[i].bondDis);

    /* translation, then a rotation of the whole ligand about the origin */
    std::vector<Point3D> moved;
    for (int i = 0; i < (int)L.size(); ++i) {
        Point3D p = L.getIthPoints(i);
        moved.push_back(Point3D(p[0] + 4.0, p[1] - 9.0, p[2] + 2.5));
    }
    L.setPoints(moved);
    const std::vector<BondInfo>& shifted = L.getBondInfoes();
    for (int i = 0; i < (int)shifted.size(); ++i) {
        CHECK_CLOSE(shifted[i].bondDis, before[i], 1e-12);
    }

    std::vector<Point3D> rotated;
    for (int i = 0; i < (int)L.size(); ++i) {
        Point3D p = L.getIthPoints(i);
        rotated.push_back(Point3D(p[0], p[2], -p[1]));   /* +90 deg about X */
    }
    L.setPoints(rotated);
    const std::vector<BondInfo>& spun = L.getBondInfoes();
    for (int i = 0; i < (int)spun.size(); ++i) {
        CHECK_CLOSE(spun[i].bondDis, before[i], 1e-12);
    }
}

/* ---- F15: ulp drift from a rebuilt rotation is not written back ---- */
TEST(f15_ulp_drift_ignored) {
    /* Flexible mode rebuilds the whole point list after each bond rotation,
     * so atoms on the rotation axis pick up ~1e-16 rounding.  Writing that
     * back to bondDis feeds an ulp into EnhancedGS and flips discrete choices
     * (measured: octane flex RMSD_lb 2.1042 -> 2.1140).  The refresh therefore
     * ignores sub-1e-9 drift but still tracks a real change. */
    std::vector<Point3D> pts;
    pts.push_back(Point3D(0.0, 0.0, 0.0));
    pts.push_back(Point3D(1.5, 0.0, 0.0));
    std::vector<std::string> types;
    types.push_back("C");
    types.push_back("C");
    std::vector<BondInfo> bonds;
    BondInfo b;
    b.aInd = 0;
    b.bInd = 1;
    b.bt = SINGLE;
    b.bondDisPow2 = 2.25;
    b.bondDis = 1.5;
    bonds.push_back(b);
    Ligand L("DRIFT", pts, types, bonds);

    /* +4 ulp on the second atom: the stored length must stay bit-identical */
    std::vector<Point3D> nudged;
    nudged.push_back(Point3D(0.0, 0.0, 0.0));
    nudged.push_back(Point3D(1.5 + 8.88e-16, 0.0, 0.0));
    L.setPoints(nudged);
    CHECK(L.getBondInfoes()[0].bondDis == 1.5);

    /* a 1e-6 change is well above the threshold and must be recorded */
    std::vector<Point3D> real;
    real.push_back(Point3D(0.0, 0.0, 0.0));
    real.push_back(Point3D(1.500001, 0.0, 0.0));
    L.setPoints(real);
    CHECK_CLOSE(L.getBondInfoes()[0].bondDis, 1.500001, 1e-12);
}
TEST(f15_index_oob_guard) {
    /* Build a one-bond ligand, then hand setPoints a shorter point list so
     * the stored bond indexes fall off the end.  The loop must skip the bond
     * (and warn on stderr) instead of reading points[1] out of range.  Run
     * this under 'make unit-assert' to make an unguarded read abort loudly. */
    std::vector<Point3D> twoPts;
    twoPts.push_back(Point3D(0.0, 0.0, 0.0));
    twoPts.push_back(Point3D(1.5, 0.0, 0.0));
    std::vector<std::string> types;
    types.push_back("C");
    types.push_back("C");
    std::vector<std::string> names;
    names.push_back("C1");
    names.push_back("C2");
    std::vector<BondInfo> bonds;
    BondInfo b;
    b.aInd = 0;
    b.bInd = 1;
    b.bt = SINGLE;
    b.bondDisPow2 = 2.25;
    b.bondDis = 1.5;
    bonds.push_back(b);
    Ligand L("TWO", twoPts, types, bonds, names);

    std::vector<Point3D> onePt;
    onePt.push_back(Point3D(5.0, 5.0, 5.0));
    L.setPoints(onePt);

    CHECK(L.size() == 1);
    CHECK(L.getBondInfoes().size() == 1);
    /* the out-of-range bond keeps its previous length; nothing was read */
    CHECK_CLOSE(L.getBondInfoes()[0].bondDis, 1.5, 0.0);
}
