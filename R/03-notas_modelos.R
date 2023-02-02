library('ExploreModelMatrix')

## Con R, usamos mucho la función model.matrix() y la sintáxis de fórmula Y ~ X1 + X2
## tal como en el siguiente ejemplo:

## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat


summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))


### ExploreModelMatrix

## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## De forma interactiva usando shiny

app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)

## Ejemplo mas complejo (2)

(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))

vd <- VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)

cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

## Ejemplo con batchs, controles mas especificos (3)

(sampleData = data.frame(
  condition = factor(rep(c("ctrl_minus", "ctrl_plus",
                           "ko_minus", "ko_plus"), 3)),
  batch = factor(rep(1:6, each = 2))))

vd <- VisualizeDesign(sampleData = sampleData,
                      designFormula = ~ 0 + batch + condition,
                      textSizeFitted = 4, lineWidthFitted = 20,
                      dropCols = "conditionko_minus")
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

## Interpreta ResponseResistant.Treatmentpre del ejercicio 2.
## Puede ser útil tomar un screenshot (captura de pantalla) y
## anotarla con líneas de colores. Si haces eso, puedes incluir la imagen en tus notas.

# A la columna pre le restamos la columna post en la tabla donde la respuesta es resistente


## ¿Por qué es clave el 0 al inicio de la fórmula en el ejercicio 3?

# Para que no se tome a batch 1 como el intercepto, y tengamos varios batchs que sirven
# como uno
