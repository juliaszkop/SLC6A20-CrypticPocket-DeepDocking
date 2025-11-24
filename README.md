# **Identification of Inhibitors of SLC6A20 (SIT1)**

This repository contains all scripts, input files, and documentation related to the computational discovery of potential inhibitors targeting **SLC6A20 (SIT1)** — a membrane transporter from the SLC6 family.
The project focuses on detecting **cryptic/allosteric pockets**, preparing a compound library, performing **virtual screening**, and building a **Deep Docking** model for large-scale prediction.

---

## **1. Selection of the Protein Structure**

SLC6A20 is experimentally available in three conformational states:

* **Outward-open**
* **Occluded**
* **Inward-open**

These conformations form a continuous transport cycle.
For our study, we selected the **outward-open** state as the starting structure because:

* It is **most accessible to extracellular ligands**, increasing the likelihood of detecting previously hidden (cryptic) or allosteric pockets.
* A **homologous SLC6 family member** in the same state contains a known allosteric pocket, supporting the choice of this conformation.

### **Structure Preprocessing**

To prepare the model for molecular dynamics:

* **Highly flexible extracellular loops**, which do not contribute to internal pocket formation, were **removed**.
* Missing regions were **remodeled using SWISS-MODEL**, using the remaining part of the structure as a template to preserve the global fold.
* The original and remodeled structures were **superimposed**, showing an RMSD of **~2.5 Å** (3754 atoms), demonstrating very good structural consistency.

<img width="1365" height="625" alt="image" src="https://github.com/user-attachments/assets/5a9db7c0-f342-4746-8284-d8aeb42df7b7" />



---

## **2. Molecular Dynamics (MD) Simulations**

To uncover potential cryptic pockets, we performed **all-atom MD simulations** of the processed structure.

### **System Setup**

* The protein was **embedded into a lipid bilayer** with a composition reflecting **eukaryotic cell membranes**.
* Based on **UniProt expression data**, SLC6A20 is predominantly found in the **kidney**, so kidney-specific membrane composition was incorporated.
* The system was simulated under **physiological conditions** (temperature, ionic strength, pressure).

<img width="1351" height="621" alt="image" src="https://github.com/user-attachments/assets/bab97cfe-cd23-447c-818a-be964f28c030" />


---

## **3. MD Trajectory Analysis**

### **Local RMSD Analysis**

Local structural fluctuations were quantified using a custom script (`rmsd.py`):

* RMSD was computed for **triplets of residues**, evaluated in **sliding time windows** (e.g., every 50 ns).
* Averaged RMSD profiles were written into separate CSV files for further analysis.

Run the script via:

```bash
python3 RMSD.py
```

This analysis revealed that the **most mobile internal regions** overlap with the **known tiagabine-binding region**, validating the approach.

<img width="1283" height="655" alt="image" src="https://github.com/user-attachments/assets/e0e82487-e74b-46e5-9767-d3654e268634" />


### **Cryptic Pocket Detection (mdPocket)**

To locate transient pockets that arise during the simulation, we applied **mdPocket**:

* mdPocket maps the **frequency and geometry** of pockets formed during MD trajectories.
* A **cryptic cavity** was detected (highlighted in red in the figure).
* Although the surrounding residues differ from those in the tiagabine binding site, the region is spatially close to a known **allosteric pocket** in a homologous SLC6 transporter.

This cryptic site was selected as the **target region for docking**.

<img width="1323" height="661" alt="image" src="https://github.com/user-attachments/assets/1d88d983-a8c4-4192-a3c3-77d06440fd35" />

### **Helical Opening Analysis (Custom COM-Based Method)**

To identify the **most open conformation** associated with the cryptic pocket, we implemented a dedicated analysis script (`helix_opening.py`). The script tracks the **center-of-mass (COM) distances** between predefined transmembrane helices and computes the **mean pairwise COM separation** for every frame of the MD trajectory. This metric captures the *global opening* of the helical bundle, enabling detection of frames where the cryptic pocket is maximally exposed.
The structure corresponding to the **maximum helix separation** is automatically exported as `most_open_cryptic_pocket.pdb` for downstream visualization and docking.

Run the script via:

```bash
python3 helix_opening.py
```
---

## **4. Compound Library Preparation (ZINC Database)**

To dock potential inhibitors, we used the **ZINC** database — a large collection of purchasable 3D-ready molecules containing more than **10.8 billion compounds**.

### **Filtering Strategy**

To obtain molecules consistent with expected inhibitor size and physicochemical properties, we filtered entries based on:

* **Heavy atom count (HAC):  (H20–H29)**
* **Partition coefficient proxy (P): P100–P600**

This filtering produced a **7 billion-compound subset** suitable for screening.

We then sampled **123,282 random molecules** using:

```bash
python3 zinc_sampler.py
```

The selected subset was converted to **PDBQT** format for docking.

---

## **5. Docking Pipeline**

### **Initial Docking with AutoDock Vina**

For preliminary screening:

* We used **AutoDock Vina** for fast ligand evaluation.
* Docking was performed into the **cryptic/allosteric pocket** identified via MD.
* The results will be used as a training set for subsequent deep learning models.

(Additional technical details about box size, exhaustiveness, and scoring can be added later.)

---

## **6. Deep Docking Model**

As a final stage, the project aims to implement a **Deep Docking (DD)** model:

* A subset of ligands is docked using Vina.
* Docking scores are used to train a neural network that predicts high-scoring ligands from the full 7B compound pool.
* This approach enables screening of **billions of molecules** at a fraction of the computational cost.

---

## **Repository Structure**

```
.
├── preprocessing/          # Structure trimming, loop removal, remodeling
├── md/                     # MD input files, system setup, logs
├── analysis/
│   ├── RMSD.py            # Local RMSD script
│   ├── mdpocket/          # Pocket detection output
├── docking/
│   ├── vina/              # Initial Vina docking
│   ├── deep_docking/      # Deep Docking model (training & inference)
├── zinc/
│   ├── zinc_sampler.py    # Random sampling of filtered compounds
└── README.md
```

---

## **Contact**

For questions, suggestions, or collaboration requests, feel free to open an issue or contact the project maintainers.

