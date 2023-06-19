#' Análisis de varianza
#'
#' Realiza un análisis de varianza para datos completamente al azar con igual número de unidades por tratamiento
#'
#' @param w (tratamiento) datos del efecto del T1 sobre la(s) variable(s)
#' @param x (tratamiento) datos del efecto del T2 sobre la(s) variable(s)
#' @param y (tratamiento) datos del efecto del T3 sobre la(s) variable(s)
#' @param z (tratamiento) datos del efecto del T4 sobre la(s) variable(s)
#' @return Una lista con la Media de los tratamientos, los Grados de libertad, la Suma de Cuadrados, el valor de F calculada y los Grados de libertad para el cálculo de F tablas
#' @export
AnVa <- function(w, x, y, z){
  #Cálculos preeliminares
  SumaW <- sum(w)
  SumaX <- sum(x)
  SumaY <- sum(y)
  SumaZ <- sum(z)

  SumaTotal <- sum(SumaW + SumaX + SumaY + SumaZ)
  ni <- 4
  N <- 16
  k <- 4

  #Cálculo de medias
  wbarra <- mean(w)
  xbarra <- mean(x)
  ybarra <- mean(y)
  zbarra <- mean(z)

  #Cálculo de Mbarra
  Mbarra <- sum(wbarra + xbarra + ybarra + zbarra)/N
  #Cálculo de Suma Cuadrados de los Tratamientos (Entre las muestras)
  SCtrat <- 1 / ni * sum(SumaW^2 + SumaX^2 + SumaY^2 + SumaZ^2) - SumaTotal^2 / N

  #Cálculo de Suma Cuadrados Total (Total)
  SCTot <- sum(w^2, x^2, y^2, z^2) - SumaTotal^2 / N

  #Cálculo de Suma Cuadrados del Error (Dentro de las muestras)
  SCErr <- SCTot - SCtrat

  #Cálculo de Grados de libertad
  Gl_EntreMuestras <- k - 1
  Gl_DentroMuestras <- N - k
  Gl_Tot <- N - 1

  #Cálculo de Cuadrados Medios
  CM_EntreMuestras <- SCtrat / Gl_EntreMuestras
  CM_DentroMuestras <- SCErr / Gl_DentroMuestras

  #Cálculo del coeficiente de variación
  CV <- sqrt(CM_DentroMuestras) / Mbarra * 100

  #Cálculo de F
  Fcalc <- CM_EntreMuestras / CM_DentroMuestras

  #Cálculo de os grados de libertad de Fcrit. para comprarar con valores de tabla de Fisher.
  #Grado de libertad en el numerador
  glnum <- (k - 1)

  #Grados de libertad en el denominador
  glden <- (N - k)

  #Lista
  dfmedias <- data.frame(wbarra, xbarra, ybarra, zbarra)
  colnames(dfmedias) <- c("Mediaw", "Mediax", "Mediay", "Mediaz")
  SumaCuad <- c(SCtrat, SCTot, SCErr)
  dfGlib <- data.frame(Gl_EntreMuestras, Gl_DentroMuestras, Gl_Tot)
  colnames(dfGlib) <- c("Tratamiento", "Error", "Total")
  Fc <- c(Fcalc)
  dfGlibFcrit <- data.frame(glnum, glden)
  colnames(dfGlibFcrit) <- c("Grados de libertad numerador", "Grados de libertad denominador")
  CoefVar <- CV
  CuadMed <- c(CM_DentroMuestras, CM_EntreMuestras)
  SCtratamiento <- SCtrat
  SCmuestras <- SCErr
  SCtotal <- SCTot

  #Tabla ANVA
  dfanva <- data.frame(SCtratamiento, SCmuestras, SCtotal, CuadMed, Fc)
  colnames(dfanva) <- c("SCTrat", "SC Error", "SC tot ", "Cuadrados Medios", "F calculada")
  rownames(dfanva) <- c("Tratamiento", "Error")
    resultado <- list("TablaAnva" = dfanva,
                    "Grados de libertad F calculada" = dfGlib,
                    "Tabla de medias" = dfmedias,
                    "CoeficienteVar" = CoefVar,
                    "Grados de libertad F critica" = dfGlibFcrit)


  return(resultado)

}

