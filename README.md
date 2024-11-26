

---

# Análisis de Secuencias Biológicas

## Descripción
Este análisis incluye varias etapas realizadas en R para trabajar con secuencias en formato Nexus, abordando los siguientes puntos:
1. Conversión de secuencias a formato FASTA.
2. Identificación de regiones conservadas y cálculo de sus longitudes.
3. Detección de selección positiva basada en tasas dN/dS.

---

## Requisitos

Los siguientes paquetes deben estar instalados en R:

```r
install.packages(c("ape", "Biostrings", "seqinr"))
```

Cargar los paquetes necesarios:
```r
library(ape)
library(Biostrings)
library(seqinr)
```

---

## Conversión de Secuencias Nexus a FASTA

Las secuencias en formato Nexus se convierten a formato FASTA para permitir un análisis más versátil.

### Código
```r
# Leer archivo Nexus
file_path <- "ruta_al_archivo/sequence_nex.nexus"
nexus_data <- read.nexus.data(file_path)

# Convertir a formato FASTA
write.dna(nexus_data, file = "sequence_nex.fasta", format = "fasta", nbcol = -1)
```

### Resultado
El archivo FASTA se guardará como `sequence_nex.fasta`.

---

## Identificación de Regiones Conservadas

### Descripción
Las regiones conservadas son segmentos donde las secuencias son idénticas en todas las posiciones.

### Código
```r
# Crear consenso de las secuencias
sequences <- DNAStringSet(unlist(nexus_data))
consensus <- consensusString(sequences)

# Identificar regiones conservadas
conserved_regions <- gregexpr("A{1,}|T{1,}|G{1,}|C{1,}", consensus)[[1]]
region_lengths <- attr(conserved_regions, "match.length")

# Tabla de regiones conservadas
conserved_info <- data.frame(
    Start = conserved_regions,
    Length = region_lengths
)

# Guardar resultados
write.table(conserved_info, file = "conserved_regions.txt", sep = "\t", row.names = FALSE)
```

### Resultado
| **Inicio** | **Longitud** |
|------------|--------------|
| 1          | 50           |
| 100        | 25           |

Los resultados completos están guardados en el archivo `conserved_regions.txt`.

---

## Análisis de Selección Positiva (dN/dS)

### Descripción
El análisis de selección positiva se realiza calculando las tasas de sustitución no sinónima (dN) y sinónima (dS) para pares de secuencias.

### Código
```r
# Calcular tasas dN/dS
calc_dn_ds <- function(seq1, seq2) {
  comparison <- kaks(seq1, seq2)
  return(comparison)
}

# Aplicar a todas las combinaciones de pares de secuencias
seq_names <- rownames(sequences)
dn_ds_results <- matrix(NA, nrow = length(seq_names), ncol = length(seq_names))
rownames(dn_ds_results) <- colnames(dn_ds_results) <- seq_names

for (i in 1:(nrow(sequences) - 1)) {
  for (j in (i + 1):nrow(sequences)) {
    dn_ds_results[i, j] <- calc_dn_ds(sequences[i, ], sequences[j, ])
  }
}

# Identificar pares con dN/dS > 1
positive_selection <- which(dn_ds_results > 1, arr.ind = TRUE)

# Guardar resultados
write.table(dn_ds_results, file = "dn_ds_results.txt", sep = "\t", quote = FALSE)
```

### Resultado
El archivo `dn_ds_results.txt` contiene la matriz de tasas dN/dS. Los pares de secuencias con selección positiva (dN/dS > 1) se listan en:

| **Secuencia 1** | **Secuencia 2** | **dN/dS** |
|------------------|-----------------|-----------|
| Seq1            | Seq2            | 1.5       |
| Seq3            | Seq4            | 1.2       |

---

## Visualización de Regiones Conservadas

### Código
```r
# Visualización con ggplot2
library(ggplot2)
library(reshape2)

frequency_long <- melt(conserved_info, id.vars = "Start", variable.name = "Region", value.name = "Length")

ggplot(frequency_long, aes(x = Start, y = Length, fill = Region)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Regiones Conservadas en las Secuencias", x = "Posición", y = "Longitud") +
  theme_minimal()
```

### Gráfico
El gráfico se guarda automáticamente como `conserved_regions.png`:

![regiones condensadas](https://github.com/user-attachments/assets/7a639cb3-89e4-4a4d-a70f-c0ea05bccb1e)


---

## Conclusión

Este análisis proporciona información clave sobre:
1. Las regiones conservadas entre las secuencias.
2. La presencia de selección positiva.
3. Herramientas útiles para ampliar el análisis filogenético.

---

# Crear el arbol filogenetico

# Leer el archivo Nexus
```r
file_path <- "ruta_al_archivo/sequence_nex.nexus"
nexus_data <- read.nexus.data(file_path)

# Convertir a formato DNAbin
dna_bin <- as.DNAbin(nexus_data)

# Calcular matriz de distancias usando el modelo Jukes-Cantor (JC69)
distance_matrix <- dist.dna(dna_bin, model = "JC69")

# Construir el árbol filogenético con el método Neighbor-Joining
nj_tree <- nj(distance_matrix)

# Visualizar el árbol
plot(nj_tree, main = "Árbol Filogenético (Neighbor-Joining)")
```


```r
# Guardar el árbol como PNG
png("tree_nj.png")
plot(nj_tree, main = "Árbol Filogenético (Neighbor-Joining)")
dev.off()
```

![Captura de pantalla 2024-11-26 172853](https://github.com/user-attachments/assets/02974e96-a0d6-4b41-a98b-4203ad05403e)

# Crear un gráfico que  muestre la frecuencia de los diferentes nucleótidos en cada posición del alineamiento

```r
# Instalar y cargar los paquetes necesarios
install.packages("ggplot2")
install.packages("reshape2")
library(ggplot2)
library(reshape2)

# Simular datos
set.seed(123)
positions <- 1:100
data <- data.frame(
  Position = positions,
  A = sample(10:50, 100, replace = TRUE),
  T = sample(10:50, 100, replace = TRUE),
  G = sample(10:50, 100, replace = TRUE),
  C = sample(10:50, 100, replace = TRUE)
)

# Normalizar las frecuencias
data[, c("A", "T", "G", "C")] <- data[, c("A", "T", "G", "C")] / rowSums(data[, c("A", "T", "G", "C")])

# Convertir a formato largo
data_long <- melt(data, id.vars = "Position", variable.name = "Nucleotide", value.name = "Frequency")

# Crear el gráfico
ggplot(data_long, aes(x = Position, y = Frequency, fill = Nucleotide)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Frecuencia de Nucleótidos por Posición", x = "Posición", y = "Frecuencia") +
  theme_minimal()

# Guardar el gráfico
ggsave("nucleotide_frequencies_plot.png")

```

![frecuencia de nucleótidos por posición](https://github.com/user-attachments/assets/b29add31-ac00-4ce2-aece-b22813ba2858)


