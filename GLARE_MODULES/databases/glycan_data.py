GLYCAN_TYPES = {

    "N":      "N-glycan",
    "O":      "O-glycan",
    "free":   "Free oligosaccharide",
    "lipid":  "Glycolipid",
    "repeat": "Repeat / polysaccharide",

}

SYSTEMS = [

    ("Heart", [
        "heart",          
        "myocardium",      
        "aorta",         
        "cardiac_myocytes", 
        "pericardium",    
    ]),

    ("Eye", [
        "eye",                
        "cornea",               
        "corneal_epithelium",   
        "retina",                
        "lens",                   
        "ciliary_body",           
        "lacrimal_gland",         
    ]),

    ("Skeletal", [
        "skeletal",        
        "cartilage",      
        "bone",              
        "synovial_joint",   
    ]),

    ("Muscle", [
        "muscle",          
        "skeletal_muscle", 
        "smooth_muscle",   
        "cardiac_muscle",  
    ]),
]

TISSUES = {

    "heart": {
        "display": "Heart (all)",
        "parent":  None,
        "terms": [
            "heart",
            "myocardium",
            "aorta",
            "cardiac_muscle_tissue",
            "pericardium",
            "regular_atrial_cardiac_myocyte",
            "regular_ventricular_cardiac_myocyte",
        ],
        "description": (
            "Entire heart and great vessels. "
            "Cardiac glycomics is dominated by N-glycans (sialylated complex-type) "
            "on membrane proteins; O-glycans are sparse. "
            "Use the sub-tissue keys for cell-type-specific queries."
        ),
    },

    "myocardium": {
        "display": "Myocardium",
        "parent":  "heart",
        "terms": [
            "myocardium",
            "cardiac_muscle_tissue",
            "heart",
        ],
        "description": (
            "Working heart muscle. Rich in sialylated N-glycans on ion channels "
            "and structural glycoproteins; O-glycans on titin and dystroglycan."
        ),
    },

    "aorta": {
        "display": "Aorta",
        "parent":  "heart",
        "terms": ["aorta"],
        "description": (
            "Large elastic artery. Proteoglycans (chondroitin/dermatan sulfate, heparan sulfate) "
            "dominate the extracellular matrix; endothelial O-glycans include Core 1 on "
            "CD34 and PECAM-1."
        ),
    },

    "cardiac_myocytes": {
        "display": "Cardiac Myocytes",
        "parent":  "heart",
        "cell_type": "Cardiomyocyte",
        "terms": [
            "regular_atrial_cardiac_myocyte",
            "regular_ventricular_cardiac_myocyte",
        ],
        "description": (
            "Isolated atrial and ventricular myocytes. "
            "N-glycan dataset from cell-surface proteomics; very few O-glycans annotated."
        ),
    },

    "pericardium": {
        "display": "Pericardium",
        "parent":  "heart",
        "terms": ["pericardium"],
        "description": (
            "Fibro-serous sac surrounding the heart. "
            "Glycan data sparse; dominated by N-glycans from fibronectin and collagen."
        ),
    },


    "eye": {
        "display": "Eye (all)",
        "parent":  None,
        "terms": [
            "cornea",
            "retina",
            "ciliary_body",
            "lens_of_camera-type_eye",
            "pigmented_layer_of_retina",
            "eye",
            "secretion_of_lacrimal_gland",
        ],
        "description": (
            "Full ocular dataset: cornea, retina, lens, ciliary body, and lacrimal secretions. "
            "Use 'corneal_epithelium' for the MUC1/MUC4/MUC16 surface and Core 1 Gal-3 site."
        ),
    },

    "cornea": {
        "display": "Cornea",
        "parent":  "eye",
        "terms": [
            "cornea",
        ],
        "description": (
            "Cornea and anterior corneal surface. Contains both stromal N-glycans "
            "(from keratocytes) and epithelial O-glycans (from surface epithelial cells). "
            "MUC1, MUC4, and MUC16 are the dominant transmembrane mucins."
        ),
    },

    "corneal_epithelium": {
        "display": "Corneal Epithelium",
        "parent":  "eye",
        "cell_type": "Non-keratinised stratified squamous epithelial cell",
        "terms": [
            "cornea",        
        ],
        "glycan_focus": "O-glycan",
        "primary_glycan": "Core 1 O-glycan — Gal(b1-3)GalNAc on MUC1 Thr86",
        "primary_lectin": "Galectin-3 (LGALS3)",
        "description": (
            "[2.1 Where] Apical non-keratinised epithelial cells of the corneal surface. "
            "We search df_glycan for all entries annotated to eye, cornea, tears, lacrimal gland, "
            "lens, and conjunctiva. This returns glycan sequences across tissue sample categories "
            "in glycowork; the 'cornea' term captures the surface epithelial dataset that includes "
            "MUC1/MUC4/MUC16 O-glycan data."
            "\n\n"
            "[2.2 What glycan type] The corneal epithelium is heavily mucin-glycosylated. "
            "MUC1, MUC4, and MUC16 are all expressed on this surface; together they dominate "
            "the glycocalyx. O-glycans constitute a substantial minority of the corneal glycome "
            "(compute exact fraction at runtime from df_glycan), with N-glycans making up the "
            "majority on stromal and membrane proteins. For lectin targeting, O-glycans are the "
            "relevant fraction because they are apically surface-exposed on mucin tandem-repeat "
            "domains."
            "\n\n"
            "[2.3 Which specific glycan] Among the O-glycans present we applied three criteria: "
            "(i) naturally abundant on the corneal epithelium, (ii) surface-exposed on the "
            "outermost mucin layer, and (iii) high binding affinity for Gal-3 in the glycan "
            "array (glycan_binding). Core 1 T-antigen (Gal(b1-3)GalNAc) on MUC1 Thr86 meets "
            "all three: it is the dominant O-glycan on MUC1 in non-cancerous corneal "
            "epithelium, it is positioned at the distal tip of the MUC1 extracellular domain "
            "where it is accessible to soluble lectins, and it is a confirmed primary ligand "
            "for Gal-3 (LGALS3) with Kd in the low-micromolar range."
            "\n\n"
            "[2.4 Disease-specific — Future Direction] In dry eye disease and ocular surface "
            "inflammation, O-glycan processing is dysregulated: sialylation decreases, "
            "truncated glycoforms (Tn antigen GalNAc, sialyl-Tn) accumulate, and Core 1 "
            "exposure increases. This enhances Gal-3 binding and drives mucin aggregation "
            "and epithelial barrier dysfunction. A therapeutic fusion protein anchored via "
            "Gal-3 to Core 1 on MUC1 could exploit this disease-state amplification. "
            "Glycan-disease signature comparison for dry eye is implemented in Section 8."
        ),
    },

    "retina": {
        "display": "Retina",
        "parent":  "eye",
        "terms": [
            "retina",
            "pigmented_layer_of_retina",
        ],
        "description": (
            "Neural retina and retinal pigment epithelium (RPE). "
            "Rich in gangliosides (GM1, GM2, GM3, GD1a) in photoreceptors; "
            "N-glycans on visual cycle proteins. MUC20 expressed in RPE."
        ),
    },

    "lens": {
        "display": "Lens",
        "parent":  "eye",
        "terms": ["lens_of_camera-type_eye"],
        "description": (
            "Crystalline lens. Predominantly N-glycans on αA-crystallin and αB-crystallin; "
            "O-glycans are minor. Unusual lipid composition (dihydrosphingomyelin)."
        ),
    },

    "ciliary_body": {
        "display": "Ciliary Body",
        "parent":  "eye",
        "terms": ["ciliary_body"],
        "description": (
            "Ciliary body smooth muscle and epithelium. "
            "Aqueous humour production; glycan data sparse in the database."
        ),
    },

    "lacrimal_gland": {
        "display": "Lacrimal Gland",
        "parent":  "eye",
        "terms": ["secretion_of_lacrimal_gland"],
        "description": (
            "Lacrimal gland secretions (tear film). "
            "Secreted O-glycoproteins including lacritin and tear lipocalin; "
            "MUC5AC and MUC7 contribute to the aqueous mucin layer."
        ),
    },

    "skeletal": {
        "display": "Skeletal (all)",
        "parent":  None,
        "terms": [
            "articular_cartilage_of_joint",
            "nasal_cartilage",
            "compact_bone_tissue",
            "acetabular_part_of_hip_bone",
            "bone_marrow",
            "bone_marrow_mast_cells",
            "bone_marrow_mast_cell",
            "synovial_fluid",
        ],
        "description": (
            "All skeletal and connective tissues. "
            "Dominated by proteoglycans and structural glycoproteins; "
            "chondroitin sulfate, keratan sulfate, and heparan sulfate are "
            "the dominant glycosaminoglycans."
        ),
    },

    "cartilage": {
        "display": "Cartilage",
        "parent":  "skeletal",
        "terms": [
            "articular_cartilage_of_joint",
            "nasal_cartilage",
        ],
        "description": (
            "Articular and nasal cartilage. "
            "Aggrecan, versican, and collagen II carry chondroitin sulfate and "
            "keratan sulfate chains; O-glycans on lubricin (PRG4) include Core 1."
        ),
    },

    "bone": {
        "display": "Bone",
        "parent":  "skeletal",
        "terms": [
            "compact_bone_tissue",
            "acetabular_part_of_hip_bone",
            "bone_marrow",
            "bone_marrow_mast_cells",
            "bone_marrow_mast_cell",
        ],
        "description": (
            "Compact bone and bone marrow. "
            "Osteocalcin and osteopontin carry N-glycans; "
            "bone marrow hematopoietic cells contribute a diverse O-glycan and "
            "lipid-glycan repertoire."
        ),
    },

    "synovial_joint": {
        "display": "Synovial Joint",
        "parent":  "skeletal",
        "terms": ["synovial_fluid"],
        "description": (
            "Synovial fluid from joint cavity. "
            "High O-glycan content: lubricin (PRG4) is heavily glycosylated "
            "with Core 1 and Core 2 O-glycans providing lubrication. "
            "MUC1 and MUC4 expressed in synoviocytes."
        ),
    },

 
    "muscle": {
        "display": "Muscle (all)",
        "parent":  None,
        "terms": [
            "skeletal_muscle_tissue",
            "muscle_of_leg",
            "muscle_structure",
            "cardiac_muscle_tissue",
            "stomach_smooth_muscle",
        ],
        "description": (
            "All muscle tissue types: skeletal, cardiac, and smooth. "
            "α-Dystroglycan is the principal glycosylated surface protein, "
            "carrying the matriglycan repeat (Xyl-Glc)n required for laminin binding. "
            "Defects cause congenital muscular dystrophies."
        ),
    },

    "skeletal_muscle": {
        "display": "Skeletal Muscle",
        "parent":  "muscle",
        "terms": [
            "skeletal_muscle_tissue",
            "muscle_of_leg",
            "muscle_structure",
        ],
        "description": (
            "Voluntary striated muscle. "
            "α-Dystroglycan O-mannosylation and matriglycan are critical for "
            "basement-membrane linkage; disruption → limb-girdle muscular dystrophies."
        ),
    },

    "smooth_muscle": {
        "display": "Smooth Muscle",
        "parent":  "muscle",
        "terms": [
            "stomach_smooth_muscle",
        ],
        "description": (
            "Involuntary smooth muscle. glycowork coverage is limited to stomach "
            "smooth muscle as the only available human term; gizzard data (avian) "
            "has been excluded. Results should be interpreted with caution — the "
            "stomach_smooth_muscle term may co-capture signal from adjacent mucosa."
        ),
    },

    "cardiac_muscle": {
        "display": "Cardiac Muscle",
        "parent":  "muscle",
        "terms": [
            "cardiac_muscle_tissue",
            "regular_atrial_cardiac_myocyte",
            "regular_ventricular_cardiac_myocyte",
        ],
        "description": (
            "Cardiac striated muscle fibres. "
            "Same database terms as myocardium; included here for motif queries "
            "from the Muscle domain perspective. "
            "N-glycans dominate; O-GlcNAcylation of myofilament proteins is extensive "
            "but not captured in df_glycan (it maps to intracellular O-GlcNAc)."
        ),
    },
}

MUCINS = {
    "MUC1": {
        "aliases": [
            "muc1", "mucin1", "mucin-1", "mucin 1", "ca15-3", "ca15.3",
            "ca 15-3", "episialin", "pem", "ema", "cd227", "peanut-reactive mucin",
        ],
        "gene":    "MUC1",
        "uniprot": "P15941",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "eye — cornea / conjunctiva / ocular surface",
            "gut epithelium — low-to-moderate on stomach, small intestine, colon",
            "liver — bile ducts",
            "pancreas — ductal epithelium",
            "breast / mammary gland",
            "lung — bronchial epithelium",
            "kidney — tubular epithelium",
        ],
        "db_terms": [
            "cornea", "eye", "colon", "colonic_mucosa",
            "intestinal_mucosa", "liver", "pancreas",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen, Galβ1-3GalNAc)",
            "Core 2 (GlcNAcβ1-6 branch on Core 1)",
            "Sialyl-Tn (NeuAcα2-6GalNAc)",
            "Tn antigen (GalNAc-α-Ser/Thr)",
            "Sialyl-Lewis X / Lewis Y",
        ],
        "pubmed": ["10880508", "15466474", "22552878", "28389562"],
        "glygen_url":    "https://www.glygen.org/protein/P15941",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/P15941",
        "notes": (
            "MUC1 VNTR contains 5 potential O-glycosylation sites per repeat. "
            "In cancer, truncated Core 1 / Tn / STn glycans expose the protein backbone. "
            "Galectin-3 binds MUC1 via β-galactoside moieties."
        ),
    },

    "MUC2": {
        "aliases": [
            "muc2", "mucin2", "mucin-2", "mucin 2", "intestinal mucin",
        ],
        "gene":    "MUC2",
        "uniprot": "Q02817",
        "type":    "Gel-forming (secreted)",
        "tissues": [
            "gut lining — colon / small intestine (dominant secreted mucin)",
            "stomach — goblet cells",
            "eye — conjunctiva (minor)",
        ],
        "db_terms": [
            "colon", "colonic_mucosa", "intestinal_mucosa",
            "large_intestine", "small_intestine", "mucosa_of_large_intestine",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan (minor)"],
        "core_motifs": [
            "Core 1 (T-antigen)",
            "Core 2",
            "Core 3 (GlcNAcβ1-3GalNAc)",
            "Sialyl-Tn",
            "Sulfated Lewis A",
        ],
        "pubmed": ["11696586", "17609482", "22012018"],
        "glygen_url":    "https://www.glygen.org/protein/Q02817",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q02817",
        "notes": "Primary structural mucin of colonic mucus gel; mouse Muc2 knockout → spontaneous colitis.",
    },

    "MUC4": {
        "aliases": [
            "muc4", "mucin4", "mucin-4", "mucin 4",
        ],
        "gene":    "MUC4",
        "uniprot": "Q99102",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "eye — conjunctiva / ocular surface",
            "gut lining — colon, small intestine",
            "lung — tracheobronchial epithelium",
            "stomach",
        ],
        "db_terms": [
            "cornea", "eye", "colon", "intestinal_mucosa",
            "mucosa_of_large_intestine",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)",
            "Core 2",
            "Sialyl-Lewis X",
        ],
        "pubmed": ["10419484", "20805372", "25677085"],
        "glygen_url":    "https://www.glygen.org/protein/Q99102",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q99102",
        "notes": "Large tandem-repeat domain (ASGP-1 is the rat name); O-glycans are densely packed on the repeat region.",
    },

    "MUC5AC": {
        "aliases": [
            "muc5ac", "mucin5ac", "mucin-5ac", "mucin 5ac", "mucin5a", "ocular muc5ac",
        ],
        "gene":    "MUC5AC",
        "uniprot": "P98088",
        "type":    "Gel-forming (secreted)",
        "tissues": [
            "eye — conjunctival goblet cells (ocular surface)",
            "gut lining — stomach (dominant gastric mucin)",
            "lung — goblet cells",
        ],
        "db_terms": [
            "eye", "cornea", "mucosa_of_stomach", "stomach_smooth_muscle",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)",
            "Core 2",
            "Sialyl-Tn",
            "Lewis B / H type 1",
        ],
        "pubmed": ["11696586", "20805372", "25677085"],
        "glygen_url":    "https://www.glygen.org/protein/P98088",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/P98088",
        "notes": (
            "Co-secreted with MUC5B in the lung; in the eye it is the primary"
            " gel mucin of the conjunctival tear film."
        ),
    },

    "MUC5B": {
        "aliases": [
            "muc5b", "mucin5b", "mucin-5b", "mucin 5b",
        ],
        "gene":    "MUC5B",
        "uniprot": "Q9HC84",
        "type":    "Gel-forming (secreted)",
        "tissues": [
            "gut lining — salivary glands / oral mucosa",
            "lung — submucosal glands",
            "cervix — endocervical mucus",
            "eye — lacrimal gland secretions (tear film)",
        ],
        "db_terms": [
            "mucosa", "mucosa_of_stomach",
            "secretion_of_lacrimal_gland",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1", "Core 2",
            "Sialyl-Lewis X", "Sialyl-Lewis A",
            "Fucosyl-Core 2",
        ],
        "pubmed": ["12499363", "17609482"],
        "glygen_url":    "https://www.glygen.org/protein/Q9HC84",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q9HC84",
        "notes": "High-molecular-weight mucin in saliva, respiratory secretions, and the tear film.",
    },

    "MUC6": {
        "aliases": [
            "muc6", "mucin6", "mucin-6", "mucin 6", "gastric mucin",
        ],
        "gene":    "MUC6",
        "uniprot": "Q9UKX1",
        "type":    "Gel-forming (secreted)",
        "tissues": [
            "gut lining — stomach pyloric antrum / duodenal Brunner's glands",
        ],
        "db_terms": [
            "mucosa_of_stomach", "duodenum",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2", "Core 3",
            "Sialyl-Tn", "Sulfated Lewis A",
        ],
        "pubmed": ["11696586", "22012018"],
        "glygen_url":    "https://www.glygen.org/protein/Q9UKX1",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q9UKX1",
        "notes": "Major mucin of pyloric antrum; marks the deep mucous neck cells.",
    },

    "MUC7": {
        "aliases": [
            "muc7", "mucin7", "mucin-7", "mucin 7", "mg2", "salivary mucin",
        ],
        "gene":    "MUC7",
        "uniprot": "Q8TAX7",
        "type":    "Non-gel-forming (secreted)",
        "tissues": [
            "gut lining — salivary glands (parotid / sublingual / submandibular)",
        ],
        "db_terms": [         
            "mucosa",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan (minor)"],
        "core_motifs": [
            "Core 1", "Core 2",
            "Sialyl-Lewis X",
        ],
        "pubmed": ["9030612", "17609482"],
        "glygen_url":    "https://www.glygen.org/protein/Q8TAX7",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q8TAX7",
        "notes": (
            "Small salivary mucin; important for oral lubrication and pathogen aggregation. "
            "Note: db_terms uses the broad 'mucosa' fallback — MUC7-specific salivary gland "
            "queries will have low tissue specificity until a dedicated glycowork term is available."
        ),
    },

    "MUC13": {
        "aliases": [
            "muc13", "mucin13", "mucin-13", "mucin 13",
        ],
        "gene":    "MUC13",
        "uniprot": "Q9H3R2",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "gut lining — colon / small intestine",
            "stomach — gastric epithelium",
            "kidney",
        ],
        "db_terms": [
            "colon", "colonic_mucosa",
            "intestinal_mucosa", "large_intestine",
        ],
        "glycan_types": ["O-glycan", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2",
        ],
        "pubmed": ["12477932"],
        "glygen_url":    "https://www.glygen.org/protein/Q9H3R2",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q9H3R2",
        "notes": "Up-regulated in colorectal cancer; may interact with ErbB2/ErbB3 signalling.",
    },

    "MUC16": {
        "aliases": [
            "muc16", "mucin16", "mucin-16", "mucin 16",
            "ca125", "ca-125", "ca 125",
            "ovarian mucin", "ocular mucin muc16",
        ],
        "gene":    "MUC16",
        "uniprot": "Q8WXI7",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "eye — cornea / conjunctiva / ocular surface (dominant)",
            "ovary — ovarian surface epithelium / peritoneum",
            "uterus / fallopian tube",
        ],
        "db_terms": [
            "cornea", "retina", "ciliary_body",
            "eye", "secretion_of_lacrimal_gland",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2",
            "Sialyl-Tn", "Tn antigen",
            "Sialyl-Lewis X", "Sialyl-Lewis A",
        ],
        "pubmed": ["15466474", "20805372", "22552878", "34379775", "30959980"],
        "glygen_url":    "https://www.glygen.org/protein/Q8WXI7",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q8WXI7",
        "notes": (
            "Serum CA-125 is a diagnostic marker for ovarian cancer. "
            "In the eye, MUC16 anchors the tear film to corneal epithelium via its SEA modules."
        ),
    },

    "MUC20": {
        "aliases": [
            "muc20", "mucin20", "mucin-20", "mucin 20",
        ],
        "gene":    "MUC20",
        "uniprot": "Q8N387",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "eye — retinal pigment epithelium",
            "kidney — proximal tubule",
            "colon",
        ],
        "db_terms": [
            "retina", "pigmented_layer_of_retina", "colon",
        ],
        "glycan_types": ["O-glycan", "N-glycan"],
        "core_motifs": [
            "Core 1", "Sialyl-Tn",
        ],
        "pubmed": ["16533819"],
        "glygen_url":    "https://www.glygen.org/protein/Q8N387",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q8N387",
        "notes": "Inhibits MET receptor phosphorylation; expressed in retinal pigment epithelium.",
    },

    "MUC3A": {
        "aliases": [
            "muc3a", "muc3", "mucin3a", "mucin-3a", "mucin 3a", "mucin3", "mucin-3",
        ],
        "gene":    "MUC3A",
        "uniprot": "Q02505",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "gut lining — small intestine / colon (enterocytes)",
        ],
        "db_terms": [
            "intestinal_mucosa", "small_intestine", "colon",
            "mucosa_of_small_intestine",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2",
            "Sialyl-Lewis X",
        ],
        "pubmed": ["10862761"],
        "glygen_url":    "https://www.glygen.org/protein/Q02505",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q02505",
        "notes": "Intestinal brush-border mucin; tethers mucus to the enterocyte apical surface.",
    },

    "MUC8": {
        "aliases": [
            "muc8", "mucin8", "mucin-8", "mucin 8",
        ],
        "gene":    "MUC8",
        "uniprot": "Q9UKN1",
        "type":    "Gel-forming (secreted)",
        "tissues": [
            "lung — trachea / bronchial submucosal glands",
            "sinus — sinonasal epithelium",
        ],
        "db_terms": [
            "lung", "trachea", "bronchial_epithelium",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2",
            "Sialyl-Tn", "Sialyl-Lewis X",
        ],
        "pubmed": ["10702290", "15955827"],
        "glygen_url":    "https://www.glygen.org/protein/Q9UKN1",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q9UKN1",
        "notes": (
            "Expressed in trachea and bronchial submucosal glands; distinct from MUC5B. "
            "Elevated in chronic sinusitis and COPD airway secretions."
        ),
    },

    "MUC12": {
        "aliases": [
            "muc12", "mucin12", "mucin-12", "mucin 12",
        ],
        "gene":    "MUC12",
        "uniprot": "Q8WXX0",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "gut lining — colon / stomach (epithelium)",
            "pancreas",
        ],
        "db_terms": [
            "colon", "colonic_mucosa", "mucosa_of_stomach",
            "large_intestine", "pancreas",
        ],
        "glycan_types": ["O-glycan", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2", "Tn antigen",
        ],
        "pubmed": ["11973330", "28389562"],
        "glygen_url":    "https://www.glygen.org/protein/Q8WXX0",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q8WXX0",
        "notes": (
            "Transmembrane mucin co-expressed with MUC2 in colonic epithelium. "
            "Emerging colorectal cancer biomarker; often upregulated alongside MUC4."
        ),
    },

    "MUC17": {
        "aliases": [
            "muc17", "mucin17", "mucin-17", "mucin 17",
        ],
        "gene":    "MUC17",
        "uniprot": "Q685J3",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "gut lining — small intestine brush border (duodenum / jejunum)",
        ],
        "db_terms": [
            "small_intestine", "duodenum", "jejunum",
            "intestinal_mucosa", "mucosa_of_small_intestine",
        ],
        "glycan_types": ["O-glycan (dominant)", "N-glycan"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Core 2",
            "Sialyl-Lewis X",
        ],
        "pubmed": ["15361234", "22819913"],
        "glygen_url":    "https://www.glygen.org/protein/Q685J3",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/Q685J3",
        "notes": (
            "Major transmembrane mucin of the intestinal brush border; functionally analogous "
            "to MUC3A but predominantly in the proximal small intestine. Protects enterocytes "
            "from luminal bacteria via its extended glycocalyx."
        ),
    },

    "MUC21": {
        "aliases": [
            "muc21", "mucin21", "mucin-21", "mucin 21",
        ],
        "gene":    "MUC21",
        "uniprot": "A6NMB1",
        "type":    "Transmembrane (type I)",
        "tissues": [
            "lung — trachea / upper airway",
            "thymus",
        ],
        "db_terms": [
            "lung", "trachea", "thymus",
        ],
        "glycan_types": ["O-glycan (dominant)"],
        "core_motifs": [
            "Core 1 (T-antigen)", "Sialylated O-glycans",
        ],
        "pubmed": ["16455651"],
        "glygen_url":    "https://www.glygen.org/protein/A6NMB1",
        "glycosmos_url": "https://glycosmos.org/glycoproteins/show/uniprot/A6NMB1",
        "notes": (
            "Membrane-tethered mucin predominantly in the upper airway epithelium and thymus. "
            "Implicated in airway defense; upregulated in lung cancer cell lines."
        ),
    },
}

TISSUE_TO_MUCINS = {

    "heart":            ["MUC1"],
    "myocardium":       ["MUC1"],
    "aorta":            [],
    "cardiac_myocytes": [],
    "pericardium":      [],

    "eye":                  ["MUC1", "MUC2", "MUC4", "MUC5AC", "MUC5B",
                             "MUC7", "MUC16", "MUC20"],
    "cornea":               ["MUC1", "MUC4", "MUC16", "MUC5AC"],
    "corneal_epithelium":   ["MUC1", "MUC4", "MUC16", "MUC5AC"],
    "retina":               ["MUC20"],
    "lens":                 [],
    "ciliary_body":         [],
    "lacrimal_gland":       ["MUC5AC", "MUC7"],

    "skeletal":         [],
    "cartilage":        [],
    "bone":             [],
    "synovial_joint":   ["MUC1", "MUC4"],

    "muscle":           [],
    "skeletal_muscle":  [],
    "smooth_muscle":    [],
    "cardiac_muscle":   ["MUC1"],
}

MOTIFS = {
    "core1": {"label": "Core 1 (T-antigen)", "pattern": r"Gal\(b1-3\)GalNAc"},
    "core2": {"label": "Core 2", "pattern": r"GlcNAc\(b1-6\)\[Gal\(b1-3\)\]GalNAc"},
    "core3": {"label": "Core 3", "pattern": r"GlcNAc\(b1-3\)GalNAc"},
    "core4": {"label": "Core 4", "pattern": r"GlcNAc\(b1-6\)\[GlcNAc\(b1-3\)\]GalNAc"},
    "tn": {"label": "Tn antigen (GalNAc-α) [heuristic]", "pattern": r"(?<!\))GalNAc(?!\()"},
    "stn": {"label": "Sialyl-Tn", "pattern": r"Neu5Ac\(a2-6\)GalNAc"},
    "lex": {"label": "Lewis X (CD15)", "pattern": r"Gal\(b1-4\)\[Fuc\(a1-3\)\]GlcNAc"},
    "lea": {"label": "Lewis A", "pattern": r"Gal\(b1-3\)\[Fuc\(a1-4\)\]GlcNAc"},
    "ley": {"label": "Lewis Y",  "pattern": r"Fuc\(a1-2\)Gal\(b1-4\)\[Fuc\(a1-3\)\]GlcNAc"},
    "slx": {"label": "Sialyl-Lewis X (sLeX)", "pattern": r"Neu5Ac.*Fuc.*GlcNAc"},
    "slea": {"label": "Sialyl-Lewis A (CA19-9)", "pattern": r"Neu5Ac\(a2-3\)Gal\(b1-3\)\[Fuc\(a1-4\)\]GlcNAc"},
    "blood_h": {"label": "Blood group H (Fucα1-2Gal)", "pattern": r"Fuc\(a1-2\)Gal"},
    "blood_a": {"label": "Blood group A (GalNAcα1-3H)", "pattern": r"GalNAc\(a1-3\).*Fuc\(a1-2\)Gal"},
    "blood_b": {"label": "Blood group B (Galα1-3H)", "pattern": r"Gal\(a1-3\).*Fuc\(a1-2\)Gal"},
    "sia_a23": {"label": "NeuAcα2-3Gal", "pattern": r"Neu5Ac\(a2-3\)Gal"},
    "sia_a26": {"label": "NeuAcα2-6", "pattern": r"Neu5Ac\(a2-6\)"},
    "polysia": {"label": "Polysialic acid", "pattern": r"Neu5Ac.*Neu5Ac"},
    "fuc": {"label": "Fucosylation (any)", "pattern": r"Fuc\("},
    "beta_gal": {"label": "β-galactoside", "pattern": r"Gal\(b"},
    "gm1": {"label": "GM1 ganglioside", "pattern": r"Gal\(b1-3\)GalNAc\(b1-4\)\[Neu5Ac\(a2-3\)\]Gal"},
    "gm2": {"label": "GM2 ganglioside", "pattern": r"GalNAc\(b1-4\)\[Neu5Ac\(a2-3\)\]Gal"},
    "gm3": {"label": "GM3 ganglioside", "pattern": r"Neu5Ac\(a2-3\)Gal\(b1-4\)Glc"},
    "sulfo_slx": {"label": "6-sulfo sialyl-Lewis X", "pattern": r"Neu5Ac\(a2-3\)Gal\(b1-4\)\[Fuc\(a1-3\)\]GlcNAc.*6S|"},
    "sulfo_gal": {"label": "3-sulfo-Gal", "pattern": r"Gal3S|3S.*Gal"},
    "gb3": {"label": "Gb3 (globotriaosylceramide)", "pattern": r"Gal\(a1-4\)Gal\(b1-4\)Glc"},
    "alpha_gal": {"label": "α-Gal (Galα1-3Gal)", "pattern": r"Gal\(a1-3\)Gal"},
}

GLYCAN_ALIASES: dict[str, str] = {

    "core1":              "Gal(b1-3)GalNAc",
    "core 1":             "Gal(b1-3)GalNAc",
    "t antigen":          "Gal(b1-3)GalNAc",
    "t-antigen":          "Gal(b1-3)GalNAc",
    "thomsen-friedenreich":"Gal(b1-3)GalNAc",
    "tf antigen":         "Gal(b1-3)GalNAc",
    "core2":              "GlcNAc(b1-6)[Gal(b1-3)]GalNAc",
    "core 2":             "GlcNAc(b1-6)[Gal(b1-3)]GalNAc",
    "core3":              "GlcNAc(b1-3)GalNAc",
    "core 3":             "GlcNAc(b1-3)GalNAc",
    "core4":              "GlcNAc(b1-6)[GlcNAc(b1-3)]GalNAc",
    "core 4":             "GlcNAc(b1-6)[GlcNAc(b1-3)]GalNAc",
    "tn":                 "GalNAc",
    "tn antigen":         "GalNAc",
    "tn-antigen":         "GalNAc",
    "stn":                "Neu5Ac(a2-6)GalNAc",
    "sialyl-tn":          "Neu5Ac(a2-6)GalNAc",
    "sialyl tn":          "Neu5Ac(a2-6)GalNAc",
    "lacnac":             "Gal(b1-4)GlcNAc",
    "lac-nac":            "Gal(b1-4)GlcNAc",
    "n-acetyllactosamine": "Gal(b1-4)GlcNAc",
    "lewis x":            "Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "lex":                "Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "sialyl lewis x":     "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "slx":                "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "slex":               "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "lewis a":            "Gal(b1-3)[Fuc(a1-4)]GlcNAc",
    "lea":                "Gal(b1-3)[Fuc(a1-4)]GlcNAc",
    "sialyl lewis a":     "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc",
    "slea":               "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc",
    "man5":               "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "man-5":              "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "high mannose":       "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
    "neu5ac":             "Neu5Ac",
    "nana":               "Neu5Ac",
    "sia a2-3":           "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
    "sia a2-6":           "Neu5Ac(a2-6)Gal(b1-4)GlcNAc",
    "alpha2-3 sialic":    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
    "alpha2-6 sialic":    "Neu5Ac(a2-6)Gal(b1-4)GlcNAc",
    "gal":                "Gal",
    "galactose":          "Gal",
    "b-galactoside":      "Gal(b1-4)GlcNAc",
    "beta-galactoside":   "Gal(b1-4)GlcNAc",
}

TARGET_CLASSES: dict[str, dict] = {

    "on_target": {
        "label": "On-target",
        "color": "#2ca02c",
        "motifs": ["Gal(b1-3)GalNAc", "Gal(b1-4)GlcNAc", "Gal(b1-3)GlcNAc",
                   "Gal(b1-4)Glc", "Lactose"],
        "description": "Core 1 / beta-galactoside -- primary Gal-3 ligands",
    },

    "off_sialyl": {
        "label": "Off-target: sialylated",
        "color": "#ff7f0e",
        "motifs": ["Neu5Ac(a2-3)Gal", "Neu5Ac(a2-6)", "Neu5Gc"],
        "description": "Sialylation caps galactoside and blocks Gal-3 binding",
    },

    "off_fucosyl": {
        "label": "Off-target: fucosylated Lewis",
        "color": "#9467bd",
        "motifs": ["Fuc(a1-2)Gal", "Fuc(a1-3)", "Fuc(a1-4)", "Lewis X", "Lewis A",
                   "Lewis Y", "Lewis B"],
        "description": "Fucosylated Lewis antigens -- structural competitors",
    },

    "off_sulfated": {
        "label": "Off-target: sulfated",
        "color": "#d62728",
        "motifs": ["6S", "3S", "4S", "sulfo", "sulfate"],
        "description": "Sulfation sterically hinders galectin CRD",
    },

    "neg_ctrl": {
        "label": "Neg. ctrl: GlcNAc-terminating",
        "color": "#888888",
        "motifs": ["GlcNAc(b1-4)GlcNAc", "GlcNAc(b1-2)Man", "GlcNAcb"],
        "description": "GlcNAc-terminating -- poor / no Gal-3 binding",
    },

}

DISEASE_SIGNATURES: dict[str, dict] = {

    "dry_eye": {
        "tissue_key": "eye",
        "label": "Dry Eye / Ocular Surface Disease",

        "signature_motifs": [
            "Neu5Ac(a2-6)GalNAc", 
            "Gal(b1-3)GalNAc"
        ],

        "output_prefix": "dry_eye",
        "refs": [
            "https://doi.org/10.1167/iovs.08-2080",
            "https://doi.org/10.1167/iovs.11-8609",
            "https://doi.org/10.1016/j.exer.2013.08.001",
        ],
    },

    "vaginal_cancer": {
        "tissue_key": "smooth_muscle",
        "label": "Vaginal / Cervical Cancer",

        "signature_motifs":  [
            "Gal(b1-4)[Fuc(a1-3)]GlcNAc",
            "GlcNAc(b1-6)[Gal(b1-3)]GalNAc",
            "Gal(b1-3)GalNAc",
            "GlcNAc(b1-4)GlcNAc",
            "Gal(b1-4)GlcNAc"
        ],

        "output_prefix": "vaginal_cancer",
        "refs": [
            "https://doi.org/10.1002/ijc.30920",
            "https://doi.org/10.1002/gcc.22744",
        ],
    },

    "kidney_cancer": {
        "tissue_key": "bone",
        "label": "Kidney / Renal Cell Carcinoma",

        "signature_motifs":  [
            "Fuc(a1-6)GlcNAc",
            "Man(a1-",
            "GlcNAc(b1-4)GlcNAc",
        ],

        "output_prefix": "kidney_cancer",
        "refs": [
            "https://doi.org/10.1021/pr500410y",
            "https://doi.org/10.1021/acs.jproteome.0c00348",
        ],
    },

    "lung_cancer": {
        "tissue_key": "smooth_muscle",
        "label": "Respiratory / Lung Cancer",

        "signature_motifs":  [
            "Fuc(a1-6)GlcNAc",
            "Gal(b1-4)[Fuc(a1-3)]GlcNAc",
            "GlcNAc(b1-4)GlcNAc",
            "GlcNAc(b1-6)[Gal(b1-3)]GalNAc"
        ],

        "output_prefix": "lung_cancer",
        "refs": [
            "https://doi.org/10.1021/pr300094t",
            "https://doi.org/10.1021/acs.jproteome.6b00693",
        ],
    },
    
}