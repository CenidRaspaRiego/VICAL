/* COMMENTS ON THE APPLICATION
        Application to determine vegetation indices in agricultural parcel
        Created by INIFAP CENID RASPA 2022 
        contacto: jimenez.sergio@inifap.gob.mx
*/
 
//=======================================================*****************************************************************************
//Drawing Tools
//=======================================================
var drawingTools = Map.drawingTools();
// Don't make imports that correspond to the drawn points.
Map.drawingTools().setLinked(false);
var debug = ui.url.get('debug', true);

//=========================================================
// Import tools, styles and images
//==========================================================
//Se importa el documento donde se tienen todas las imagenes 
var imp = require('users/InifapCenidRaspa/VICAL:Exportaciones');
var imp2 = require('users/InifapCenidRaspa/VICAL:VegetationIndex');
var St = require('users/InifapCenidRaspa/VICAL:Style');
//LOGOS
//var logo_INIFAP=ee.Image('users/InifapCenidRaspa/Cohen/LogoM').resample('bicubic');

// ===============================================================================****************************************************
// 1. COLECCIÓN DE IMAGENES LANDSAT Y SENTINEL
// ===============================================================================
//Pide al servidor para estimar elevación
var elevation = ee.Image('CGIAR/SRTM90_V4');
// Pide al servidor conjunto de datos Landsat 7 y 8
function ColeccionImagenSR(fecha, recorte, umbral)
{
  var L9sr = imp.ColeccionLandsatSR(fecha, 'LC09', recorte, umbral);
  var L7sr = imp.ColeccionLandsatSR(fecha, 'LE07', recorte, umbral);
  var L8sr = imp.ColeccionLandsatSR(fecha, 'LC08', recorte, umbral);
  // Armoniza los sensores a los del sensor OLI
  var L7a = L7sr.map(imp.TMaOLI);
  // Selecciona las mismas bandas que los otros conjuntos de datos
  var L8a = L8sr.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  var L9a = L9sr.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  // Une las tres series de imágenes
  var serieT =L7a.merge(L8a).merge(L9a).sort('system:time_start');
	return serieT;
}

function ColeccionImagenSR_4(fecha, recorte, umbral)
{
  var L5sr = imp.ColeccionLandsatSR(fecha, 'LT05', recorte, umbral);
  var L4sr = imp.ColeccionLandsatSR(fecha, 'LT04', recorte, umbral);
  var L5a = L5sr.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  var L4a = L4sr.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  var serieT =L4a.merge(L5a).sort('system:time_start');
	return serieT;
}

function maskEdges(s2_img) {
  return s2_img.updateMask(
      s2_img.select('redE4').mask().updateMask(s2_img.select('B9').mask()));
}
function maskClouds(img) {
  var clouds = ee.Image(img.get('cloud_mask')).select('probability');
  var isNotCloud = clouds.lt(70);
  return img.updateMask(isNotCloud);
}

// Pide al servidor conjunto de datos sentinel
function ColeccionImagenSentinelSR(fecha, recorte, umbral)
{
  var S2sr = imp.ColeccionImagenSentinelSR(fecha, recorte, umbral);
return S2sr;
}

function ColeccionImagenAMBOS(fecha, recorte, umbral)
{
  var L8Conjunto=ColeccionImagenSR(fecha, recorte, umbral)
  var S2sr = imp.ColeccionImagenSentinelSR(fecha, recorte, umbral);
  var S2a = S2sr.map(imp.MSIaOLI);
  var serieT = S2a.merge(L8Conjunto).sort('system:time_start');
	return serieT;
}

// ===============================================================================***************************************************
// Function to download images and assign names
// ===============================================================================
function getName(imagen){
      var time_start=imagen.get('system:time_start');
      var date=ee.Date(time_start);
      var year=ee.Number(date.get('year'));
      var month=ee.Number(date.get('month'));
      var day=ee.Number(date.get('day'));
      var sensor=imagen.get('sensor');
      var satelite=imagen.get('satelite');
      var Name_image=ee.String(year).cat(ee.String('_')).cat(ee.String(month)).cat(ee.String('_')).cat(ee.String(day)).cat(ee.String('_')).cat(ee.String(satelite)).cat(ee.String('_')).cat(ee.String(sensor));
      return  Name_image;
 }
 function getURLNDVI (image, indice, recorte, escala){ 
    var Name=getName(image);
    var recorte2=image.clip(recorte.geometry());
    var NameMap=Name.cat(ee.String('_'+indice));
    var NameDownload=NameMap.getInfo();
    var downloadUrl = ((recorte2.select(indice)).multiply(10000).int()).divide(10000)
                    .getDownloadURL({ 
                      name: NameDownload,
                      scale: escala});
    return downloadUrl;
}
function getVectorFile (recorte){ 
// Export an ee.FeatureCollection as an Earth Engine asset.
    var srr=ee.FeatureCollection(recorte);
    var downloadUrl = srr
                    .getDownloadURL({ 
                      format:'kml',
                      filename: 'Polygons'});
    return downloadUrl;
}
// ===============================================================================***************************************************
// Function for other indices that require adjustment
// ===============================================================================
// Atmospherically Resistant Vegetation Index (ARVI)
function ARVI (imagen)
{
  var arvi = imagen.expression(
    'float(NIR-(R-Y*(B-R)))/(NIR+(R-Y*(B-R)))',{
    NIR : imagen.select('nir'),
    R : imagen.select('red'),
    B: imagen.select('blue'),
    Y  : parseFloat(Eps_in.getValue()),
    });
    arvi = arvi.rename(['ARVI']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile", ]);
    return arvi;
}

function RegresionF (imagen)
{
  var arvi = imagen.expression(
    'a*VI+b',{
    VI : imagen.select(nombreInd),
    a  : parseFloat(Ra_in.getValue()),
    b  : parseFloat(Rb_in.getValue()),
    });
    arvi = arvi.rename(['linealModel']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile", ]);
    return imagen.addBands(arvi);
}
function RegresionFc (imagen)
{
  var arvi = imagen.expression(
    'a+b*VI+c*VI**2',{
    VI : imagen.select(nombreInd),
    a  : parseFloat(Rac_in.getValue()),
    b  : parseFloat(Rbc_in.getValue()),
    c  : parseFloat(Rcc_in.getValue()),
    });
    arvi = arvi.rename(['QuadraticModel']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile", ]);
    return imagen.addBands(arvi);
}
function RegresionFp (imagen)
{
  var arvi = imagen.expression(
    'a*VI**float(b)',{
    VI : imagen.select(nombreInd),
    a  : parseFloat(Ra_in.getValue()),
    b  : parseFloat(Rb_in.getValue()),
    });
    arvi = arvi.rename(['PotentialModel']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile", ]);
    return imagen.addBands(arvi);
}
function RegresionFe (imagen)
{
  var arvi = imagen.expression(
    'a*exp(float(b*VI))',{
    VI : imagen.select(nombreInd),
    a  : parseFloat(Ra_in.getValue()),
    b  : parseFloat(Rb_in.getValue()),
    
    });
    arvi = arvi.rename(['ExponentialModel']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile", ]);
    return imagen.addBands(arvi);
}
function EVI (imagen)
{
  var evi = imagen.expression(
    '2.5*float((NIR-RED)/(NIR+C1*RED-C2*BLUE+L))',{
    NIR : imagen.select('nir'),
    RED : imagen.select('red'),
    BLUE: imagen.select('blue'),
    C1  : parseFloat(C1_in.getValue()),
    C2  : parseFloat(C2_in.getValue()),
    L   : parseFloat(L_in.getValue())
    });
    evi = evi.rename(['EVI']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile"]);
    return evi;
}
function EVI2 (imagen)
{
  var evi2 = imagen.expression(
    '2.5*(NIR-RED)/(NIR+C1*RED+1)',{
    NIR : imagen.select('nir'),
    RED : imagen.select('red'),
    C1:parseFloat(C12_in.getValue())
    });
    evi2 = evi2.rename(['EVI2']).float().copyProperties(imagen, 
      ["system:time_start", "satelite", "sensor", "tile"]);
    return evi2;
}
function OSAVI(imagen)
{
  var osavi = imagen.expression(
  '(1.16*(NIR - RED) / (NIR + RED + X))',
  {
    RED : imagen.select('red'),
    NIR : imagen.select('nir'),
    X   : parseFloat(X_in.getValue())
  });
  osavi = osavi.rename(['OSAVI']).copyProperties(imagen, 
    ["system:time_start", "satelite", "sensor", "tile"]);
  return osavi;
}
function SAVI (imagen)
{
  var savi = imagen.expression(
  '((1+L) * (NIR - RED) / (L + NIR + RED))',
  {
    RED : imagen.select('red'),
    NIR : imagen.select('nir'),
    L   : parseFloat(LS_in.getValue())
  });
  savi = savi.rename(['SAVI']).copyProperties(imagen, 
    ["system:time_start", "satelite", "sensor", "tile"]);
  return savi;
}
function ATSAVI (imagen)
{
  var atsavi = imagen.expression(
  '(a*(NIR-a*R-b) / (R+a*NIR-a*b+X*(1+a**2)))',
  {
    R   : imagen.select('red'),
    NIR : imagen.select('nir'),
    a   : parseFloat(a_in.getValue()),
    b   : parseFloat(b_in.getValue()),  
    X   : parseFloat(XS_in.getValue())
  });
  atsavi = atsavi.rename(['ATSAVI']).copyProperties(imagen, 
    ["system:time_start", "satelite", "sensor", "tile"]);
  return atsavi;
}
function TSAVI (imagen)
{
  var tsavi = imagen.expression(
  '(a*(NIR-a*R-b) / (R+a*NIR-a*b))',
  {
    R   : imagen.select('red'),
    NIR : imagen.select('nir'),
    a   : parseFloat(aT_in.getValue()),
    b   : parseFloat(bT_in.getValue()),  
  });
  tsavi = tsavi.rename(['TSAVI']).copyProperties(imagen, 
    ["system:time_start", "satelite", "sensor", "tile"]);
  return tsavi;
}
function WDRVI (imagen)
{
  var wdrvi = imagen.expression(
  '((a*NIR - RED) / (a*NIR + RED))',
  {
    RED : imagen.select('red'),
    NIR : imagen.select('nir'),
    a   : parseFloat(alpha_in.getValue())
  });
  wdrvi = wdrvi.rename(['WDRVI']).copyProperties(imagen, 
    ["system:time_start", "satelite", "sensor", "tile"]);
  return wdrvi;
}
function ProdRelativa(indice, poligonos)
{
  var media = indice.reduceRegions({
    collection  : poligonos,
    reducer : ee.Reducer.mean(),
    scale : 10
  });

  var mediaI = media.filter(ee.Filter.notNull(['mean']))
    .reduceToImage({
      properties: ['mean'],
      reducer: ee.Reducer.first()
  });

  indice = indice.clip(poligonos);
  var prod = indice.divide(mediaI).rename('ProdMedia');
  prod = prod.copyProperties(indice, indice.propertyNames());
  var geom = ee.Geometry(indice.get('system:footprint'));
  prod = prod.set('system:footprint', geom);
  return prod;
}

function MediaS(indice, poligonos)
{
  var media = indice.reduceRegions({
    collection  : poligonos,
    reducer : ee.Reducer.mean(),
    scale : 10
  });

  var mediaI = media.filter(ee.Filter.notNull(['mean']))
    .reduceToImage({
      properties: ['mean'],
      reducer: ee.Reducer.first()
  });

  indice = indice.clip(poligonos);
  var prod = mediaI.rename('ProdMedia');
  prod = prod.copyProperties(indice, indice.propertyNames());
  var geom = ee.Geometry(indice.get('system:footprint'));
  prod = prod.set('system:footprint', geom);
  return prod;
}

// ===============================================================================***************************************************
// 3. Funciones de Interfaz del usuario
// ===============================================================================
// Define parámetros de visualización
function Leyenda(titulo, vis)
{
  // Crea una imagen en miniatura de la barra de colores para usar 
  // en la leyenda de la paleta de colores dada
  function ParamBarraColor(paleta) {
    return {
      bbox: [0, 0, 1, 0.1],
      dimensions: '100x10',
      format: 'png',
      min: 0,
      max: 1,
      palette: paleta,
    };
  }

  // Crea la barra de color para la leyenda.
  var BarraColor = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: ParamBarraColor(vis.palette),
    style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},
  });

  // Crea un panel con tres números para la leyenda..
  var LeyendaEtiquetas = ui.Panel({
    widgets: [
      ui.Label(Math.round(100*vis.min)/100, {margin: '4px 8px'}),
      ui.Label(
          Math.round(100*(vis.max +  vis.min) / 2) / 100,
          {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
      ui.Label(Math.round(100*vis.max)/100, {margin: '4px 8px'})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });

  var LeyendaTitulo = ui.Label({
    value: titulo,
    style: {fontWeight: 'bold'}
  });

  // Retorna el panel de leyenda.
  var PanelLeyenda = ui.Panel({style: {position: 'bottom-left'}});
  PanelLeyenda.widgets().set(0, LeyendaTitulo);
  PanelLeyenda.widgets().set(1, BarraColor);
  PanelLeyenda.widgets().set(2, LeyendaEtiquetas);
  return(PanelLeyenda);
}

function percentiles(imagen, geom)
{
  var VminMax = imagen.reduceRegion(
    {
    reducer: ee.Reducer.percentile([1, 99]),
    geometry: geom,
    scale: 30,
    maxPixels: 1e9
    }); 
  
  var nombre = imagen.bandNames();
  var p1 = nombre.getString(0).cat(ee.String('_p1'));
  var p99 = nombre.getString(0).cat(ee.String('_p99'));
  var min = VminMax.getNumber(p1);
  var max = VminMax.getNumber(p99);
  
  var Rango = ee.List([min, max]);
  return Rango;
}

var sensor;
var bordes;
var punto;

// ===============================================================================****************************************************
// 4. PANEL DE NAVEGACION-TEXTO DE PRESENTACIÓN
// ===============================================================================
// Crea panel de controles
var panel = ui.Panel({style: {width: '360px'}});

// Crea panel de gráfica
var pangr = ui.Panel({style: {position: 'bottom-right', width: '400px'}});
//var logo_EE=ui.Thumbnail({image:logo_INIFAP,params:{bands:['b1','b2','b3'],min:0,max:255}, style:{width:'80%', padding: '10px 10px 10px 10px',margin: 'auto'}});
//panel.add(logo_EE);
//Agrega texto de version y titulo de APP
panel.add( ui.Label({value: 'Version 1.0 - last update: april 2022', style: St.version}));
panel.add( ui.Label({value: 'VICAL: VEGETATION INDICES CALCULATOR', style: St.titleAPP}));
// Agrega dos botones de como usar y acerca de 
var ButtonFAQ = ui.Button({label: 'How to use this tool?', style: St. ButtonPrincipal1, onClick: function(){Map.add(FAQ_PANEL)}});
var ButtonAbout = ui.Button({label: 'About VICAL', style: St.ButtonAbout,  onClick: function(){Map.add(ABOUT_PANEL)}});
var Panels_Button= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   widgets: [ButtonFAQ,ButtonAbout]
});

var VIHelp_Button = ui.Button({
  label: '?',
  style: St.HelpGeneral,
  //imageUrl:'https://fonts.google.com/icons?selected=Material%20Symbols%20Outlined%3Afile_download%3AFILL%400%3Bwght%40400%3BGRAD%400%3Bopsz%4048',
  onClick: function(){Map.add(INDICE_PANEL)}
});
var SATHelp_Button = ui.Button({
  label: '?',
  style: St.HelpGeneral,
  onClick: function(){Map.add(SAT_PANEL)}
});
var PRHelp_Button = ui.Button({
  label: '?',
  style: St.HelpGeneral,
  onClick: function(){Map.add(PR_PANEL)}
});
// ===============================================================================
// DIFERENTES CUADROS DE DIALOGOS
// ===============================================================================
var FAQ_PANEL=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
        ui.Label({ value:'¿How to use this Tool?',style: St.titleWidwet}), 
        ui.Label({ value:'0 – Locate your parcel on the map and start drawing polygons with the drawing tools.', style: St.conte}),
        ui.Label({ value:'1 – Enter a start and end date for the time series, using the format Year-month-day (xxxx-xx-xx)', style: St.conte}),
        ui.Label({ value:'2 – Select a set of images: 1) Landsat (7,8 and 9), 2) Sentinel-2, 3) LandSat (7, 8 and 9) and sentinel-2, 4) Landsat (4 and 5)', style: St.conte}),
        ui.Label({ value:'3 – Select a vegetation index (set index values if necessary)', style: St.conte}),
        ui.Label({ value:'4 – Optional: Calculate the "Weighting Factor" or map regression models (lineal, quadratic, potential and exponential model)', style: St.conte}),        
        ui.Label({ value:'5 – Click on "CALCULATE" and the vegetation index selected for the polygon appears on the map, the first image found in the date range is displayed.', style: St.conte}),
        ui.Label({ value:'6 – Different layers will be shown on the map: RGB, Vegetation Indices, Weighting Factor, regression model and polygon',style: St.conte}),
        ui.Label({ value:'7 – Click inside a polygon to plot the vegetation index time series. ',style: St.conte}),
        ui.Label({ value:'8 – Click on "New Polygon" to calculate the vegetation index of other polygons',style: St.conte}),
        ui.Label({ value:'9 – click on "Edit Polygon" to edit the geometry of the drawn polygons',style: St.conte})
    ],
    style: {position: 'top-center',  shown: true,  width: '40%', height: '50%',  padding: '5px', margin: '10px'}
  });
//Cuadro de dialogo de acerca de
var ABOUT_PANEL=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
        ui.Label({value:'About VICAL',style: St.titleWidwet}),
        ui.Label({value:'VICAL is a tool to calculate 23 vegetation indices (commonly used in agricultural applications with remote sensing) with Sentinel-2 and LandSat images of any agricultural area ',style: St.conte}),
        ui.Label({ value:'Authors:',style: {fontWeight: 'bold', textAlign: 'center',padding: '5px',  margin: 'auto'}}),
        ui.Label({ value:'INIFAP CENID-RASPA',  style: {fontWeight: 'bold', color:'#949494', textAlign: 'left',padding: '5px',  margin: 'auto'}}),
        ui.Label({ value:'Sergio Iván Jiménez-Jiménez; Mariana de Jesús Marcial-Pablo; Waldo Ojeda-Bustamante; Ernesto Sifuentes-Ibarra; Marco Antonio Inzunza-Ibarra, Ignacio Sánchez-Cohen*' , style: {fontWeight: 'normal', color:'#949494', textAlign: 'center',padding: '5px',  margin: 'auto'}}),
    ],
    style: { position: 'top-center',shown: true,width: '30%',height: '50%',padding: '0%',margin: '0%',}});  
// CUADRO DE DIALOGO DE SATELITE
var SAT_PANEL=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
        ui.Label({ value:'Satellites',style: St.titleWidwet}), 
        ui.Label({ value:'Vegetation indices can be extracted from: (1) Landsat (7, 8 and 9), (2) Sentinel-2 and (3) Landsat (7, 8 and 9) and Sentinel-2, 4) Landsat (4 and 5). VICAL used the atmospherically corrected land surface reflectance images', style: St.conte}),
        ui.Label({ value:'Sensor | Spatial Resolution (m/pixel) | Repeat Coverage (days) | Dataset Availability in GEE ', style: St.conte}),        
        ui.Label({ value:'LandSat 4 TM | 30 - 90 |  16 | 22/08/1982-24/06/1993 |', style: St.conte}),
        ui.Label({ value:'LandSat 5 TM | 30 - 90 |  16 | 16/03/1993 - 05/05/2012 |', style: St.conte}),
        ui.Label({ value:'LandSat 7 ETM+ | 30 - 60 |  16 | 01/01/1999-present |', style: St.conte}),
        ui.Label({ value:'LandSat 8 OLI| 30 - 60 | 16 | 11/04/2013- present |', style: St.conte}),
        ui.Label({ value:'LandSat 9 OLI-2| 30 - 60 | 16 | 31/10/2021- present |', style: St.conte}),
        ui.Label({ value:'Sentinel-2 MSI | 10 - 60 | 5 |  28/03/2017-present |', style: St.conte}),
        ui.Label({ value:'References ', style: {fontSize: '14px',  textAlign: 'center', margin: 'auto', fontWeight: 'bold'}}),
        ui.Label({ value:'ETM + data are spectrally adjusted to OLI spectral bands using Roy et al., (2016). ' , style: St.conte}),
        ui.Label({ value:'Roy, D. P., Kovalskyy, V., Zhang, H. K., Vermote, E. F., Yan, L., Kumar, S. S., & Egorov, A. (2016).'+ 
                          'Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation '+
                          'index continuity. Remote sensing of Environment, 185, 57-70',
           style: {fontSize:'11px', textAlign: 'left',padding: '7px',  margin: 'auto'},
           targetUrl:'http://dx.doi.org/10.1016/j.rse.2015.12.024'
        }),
        ui.Label({ value:'MSI data are spectrally adjusted to OLI spectral bands using Claverie et al., (2018). ' , style: St.conte}),
        ui.Label({ value:'Claverie, M., Ju, J., Masek, J.G., Dungan, J.L., Vermote, E.F., Roger, J.C., Skakun, S. v., Justice, C., 2018. The Harmonized Landsat and Sentinel-2 surface reflectance data set. Remote Sensing of Environment 219, 145–161.',
           style: {fontSize:'11px', textAlign: 'left',padding: '7px',  margin: 'auto'},
           targetUrl:'https://doi.org/10.1016/J.RSE.2018.09.002'
        }),
    ],
    style: {position: 'top-center',  shown: true,  width: '40%', height: '60%',  padding: '5px', margin: '10px'}
  });
var PR_PANEL=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
        ui.Label({ value:'Weighting factor (WF)',style: St.titleWidwet}), 
        ui.Label({ value:'Weighting factor (WF) is the ratio between the value of the index in a pixel and the average of the index in the polygon (parcel)', style: St.conte}),
        ui.Label({ value:'The WF of an agricultural parcel is a normalized indicator of the productive potential of each pixel in an image', style: St.conte}),        
        ui.Label({ value:'WF=(IVpixel)/IVavg', style: St.conte}),
        ui.Label({ value:'where: WF: Weighting Factor (adimensional); IVpixel: index value in a pixel, IVavg: average of the index in the polygon(parcel)', style: St.conte}),
    ],
    style: {position: 'top-center',  shown: true,  width: '30%', height: '40%',  padding: '5px', margin: '10px'}
  });
//Cuadro de Dialogo de Etapas de Ceval  
var INDICE_PANEL=ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    widgets:[
        ui.Label({ value:'VEGETATION INDICES AVAILABLE IN VICAL',
                  style: { fontSize: '14px', fontWeight: 'bold', textAlign: 'center',padding: '20px',  margin: 'auto'}
        }),
        ui.Label({ value:'The following list shows the vegetation indices and their formula ', style: {fontSize: '13px'}}),
        ui.Label({ value:' ---------------------------------------------------------------------------', style: {fontSize: '13px'}}),   
        ui.Label({ value:'1. ARVI*:  Atmospherically Resistant Vegetation Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR - RB)/(NIR + RB) where RB=R-γ(B-R); => default values: γ=1', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'2. ATSAVI*:  Type Soil Atmospheric Impedance Vegetation Index', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    [a(NIR - aR - b)]/[(R + aNIR - ab + X(1+a^2))]; => default values: a=1; b=0; X=0.08', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'3. DVI: Difference vegetation index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR-R) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'4. EVI*: Enhanced Vegetation index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    2.5((NIR - R)/(NIR + C1*R - C2*B + L)); => default values: C1=6.0; C2=7.5; L=1.0', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'5. EVI2*: Enhanced Vegetation index2 ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    2.5 ((NIR - R)/(NIR + C1*R + 1)); => default values: C1=2.4  ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'6. GNDV: Green Normalized Difference Vegetation Index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR-G)/(NIR+G) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'7. MSAVI2: Modified Soil Adjusted Vegetation Index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    ((2NIR+1)-√((2NIR+1)^2-8(NIR-R) ))/2', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'8. MSI: Moisture Stress Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    SWIR1/NIR ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'9. MTVI: Modiﬁed triangular vegetation index   ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    1.2[1.2(NIR-G)-2.5(R-G)]   ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'10. MTVI2: Modiﬁed Triangular Vegetation Index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    1.5[1.2(NIR-G)-2.5(R-G)]/√((2NIR+1)^2-(6NIR-5√R)-0.5)  ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'11. NDTI: Normalized Difference Tillage Index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (SWIR1 - SWIR2)/(SWIR1 + SWIR2) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'12. NDVI: Normalized Difference Vegetation Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR-R)/(NIR+R) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'13. NDWI: Normalized Difference Water Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (G - NIR)/(G + NIR) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'14. OSAVI*: Optimized Soil Adjusted Vegetation Index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR-R)/(NIR+R+X); => default values: X=0.16 ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'15. RDVI: Renormalized Difference Vegetation Index  ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR-R)/sqr(NIR+R) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'16. RI: Redness Index', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (R-G)/(R+G) ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'17. RVI: Ratio Vegetation Index)', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    R/NIR ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'18. SAVI*: Soil Adjusted Vegetation Index', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (NIR-R)/(NIR+R+L)(1+L); => default values: L=0.5', style: {fontSize: '13px', color:'#727171'}}),        
        ui.Label({ value:'19. TVI: Triangular Vegetation Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    0.5[120(NIR-G)-200(R-G)] ', style: {fontSize: '13px', color:'#727171'}}),
        ui.Label({ value:'20. TSAVI*: Transformed Soil Adjusted Vegetation Index', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    [a(NIR-aR-b)]/[(R+aNIR-ab)]; => default values: a=1; b=0 ', style: {fontSize: '13px', color:'#727171'}}),            
        ui.Label({ value:'21. VARI: Visible Atmospherically Resistant Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (G - R)/(G + R - B) ', style: {fontSize: '13px', color:'#727171'}}),    
        ui.Label({ value:'22. VIN: Vegetation Index Number ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    NIR/R ', style: {fontSize: '13px', color:'#727171'}}),    
        ui.Label({ value:'23. WDRVI*: Wide Dynamic Range Vegetation Index ', style: {fontSize: '13px'}}),
        ui.Label({ value:'->    (∝NIR-R)/(∝NIR+R); => default values: ∝=0.2', style: {fontSize: '13px', color:'#727171'}}),   
        ui.Label({ value:' ---------------------------------------------------------------------------', style: {fontSize: '13px'}}),  
        ui.Label({ value:'*Some coefficients can be configured.', style: {fontSize: '12px', color:'#de630e'}}),   
    ],
    style: { position: 'top-center',shown: true,width: '40%',height: '70%',padding: '0%',margin: '0%',}});  

var CloseButton = ui.Button({label: 'Close',style: St.closeButton,
  onClick: function(){
    Map.remove(FAQ_PANEL);
  }
});
var CloseButtonABOUT = ui.Button({label: 'Close',style: St.closeButton,
  onClick: function(){
    Map.remove(ABOUT_PANEL);
  }
});
var CloseButtonSAT = ui.Button({label: 'Close',style: St.closeButton,
  onClick: function(){
    Map.remove(SAT_PANEL);
  }
});
var CloseButtonINDICE = ui.Button({label: 'Close',style: St.closeButton,
  onClick: function(){
    Map.remove(INDICE_PANEL);
  }
});
var CloseButtonPR = ui.Button({label: 'Close',style: St.closeButton,
  onClick: function(){
    Map.remove(PR_PANEL);
  }
});
var contact = ui.Label({
  value: 'Contact - jimenez.sergio@inifap.gob.mx',
  style: {fontSize: '12px', textAlign: 'center',padding: '5px 5px', margin: 'auto'}
});

panel.add(Panels_Button);
panel.add(contact);

//SE AGREGAN LOS BOTONES DE CERRAR A LOS PANELES
FAQ_PANEL.add(CloseButton);
ABOUT_PANEL.add(CloseButtonABOUT);
INDICE_PANEL.add(CloseButtonINDICE);
SAT_PANEL.add(CloseButtonSAT);
PR_PANEL.add(CloseButtonPR);

// ===============================================================================************************************************
// 5. PANEL DE NAVEGACION-SELECCION DE OPCIONES
// ===============================================================================
// 5.1. Agrega panel con dos etiquetas para la fecha inicial y final y los cuadros de texto
var Box_Data= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: St.styleDatWiget,//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Start and End date', style: St.styleWiget}),]
});

var fechahoy=ee.String(ee.Date(Date.now()).get('year')).cat('-').cat(ee.Date(Date.now()).format('MM')).cat('-').cat(ee.Date(Date.now()).format('dd'));
var fechaayer=ee.String(ee.Date(Date.now()).advance(-1,'year').get('year')).cat('-').cat(ee.Date(Date.now()).format('MM')).cat('-').cat(ee.Date(Date.now()).format('dd'));

var fechain = ui.Textbox({
  placeholder: 'YYYY-MM-DD', 
  value: fechaayer.getInfo(),
  style: {maxWidth: '120px'}});
var fechafin = ui.Textbox({
  placeholder: 'YYYY-MM-DD', 
  value: fechahoy.getInfo(),
  style: {maxWidth: '120px'}});
  
panel.add(Box_Data);
panel.add(ui.Panel([fechain,fechafin], ui.Panel.Layout.flow('horizontal')));

//5.2.1 para porcentaje de nubes
var Nubes_text= ui.Label({value: 'Cloudy Pixel Percentage:',
    style: St.styleTexTConf
    });
var kh_in= ui.Textbox({
      placeholder: '30', 
      value: '30',
      style: {maxWidth: '70px'}
    });

var value_ConfHS= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   widgets: [Nubes_text,kh_in]
});
panel.add(value_ConfHS);

// 5.2. Agrega selector de índices
var indiceT= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: St.styleDatWiget,//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Select satellite and Vegetation Index', style: St.styleWiget}),]
});

var minPalete;
var maxPalete;
var Sensor = ui.Select({ 
 items: [ 
        { label:"1. LANDSAT (7, 8 and 9)", value: 0 }, 
        { label:"2. SENTINEL-2", value: 1}, 
        { label:"3. LANDSAT (7 and 8) and SENTINEL-2", value: 2}, 
        { label:"4. LANDSAT 4 and 5", value: 3}, 
 ],
  placeholder : 'Selecciona un satelite',
  value: 0,
  style: {width: '270px', border: '1px solid darkgray'},
});  
var indice = ui.Select({ 
 items: [ 
        { label:"1. ARVI",    value: 0 }, 
        { label:"2. ATSAVI",  value: 1 },         
        { label:"3. DVI",     value: 2 },    
        { label:"4. EVI",     value: 3 }, 
        { label:"5. EVI2",    value: 4 },
        { label:"6. GNDVI",   value: 5 },
        { label:"7. MSAVI2",  value: 6 },
        { label:"8. MSI",     value: 7 },
        { label:"9. MTVI",    value: 8 },
        { label:"10. MTVI2",  value: 9 },
        { label:"11. NDTI",   value: 10 },
        { label:"12. NDVI",   value: 11 },
        { label:"13. NDWI",   value: 12 },
        { label:"14. OSAVI",  value: 13 },
        { label:"15. RDVI",  value: 14 },        
        { label:"16. RI",     value: 15 },
        { label:"17. RVI",    value: 16 },
        { label:"18. SAVI",   value: 17 },
        { label:"19. TVI",    value: 18 },        
        { label:"20. TSAVI",  value: 19 },
        { label:"21. VARI",   value: 20 },  
        { label:"22. VIN",    value: 22 },
        { label:"23. WDRVI",  value: 23 },
 ],
  placeholder : 'Selecciona un indice de vegetación',
  value: 11,
  style: {width: '270px', border: '1px solid darkgray'},
  onChange: function(value) {
    if (value === 0) { //es para el ARVI
        BoxC_ARVI.style().set('shown', true);
        value_ARVI.style().set('shown', true);
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);        
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }
    else if  (value === 1) { //es para el ATSAVI
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false);        
        BoxC_ATSAVI.style().set('shown', true);
        value_ATSAVI.style().set('shown', true);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }
    else if (value==3) {//es para el EVI
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', true);
        value_EVI.style().set('shown', true);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }
    else if (value==4) {//es para el EVI2
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', true);
        value_EVI2.style().set('shown', true);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }
    else if (value==13) {//es para el OSAVI
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', true);
        value_OSAVI.style().set('shown', true);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }    
    else if (value==17) {//es para el SAVI
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', true);
        value_SAVI.style().set('shown', true);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }  
    else if (value==19) {//es para el TSAVI
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', true);
        value_TSAVI.style().set('shown', true);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }  
    else if (value==23) {//es para el wdrvi
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', true);
        value_WDRVI.style().set('shown', true);
    }  
    else {
        BoxC_ARVI.style().set('shown', false);
        value_ARVI.style().set('shown', false); 
        BoxC_ATSAVI.style().set('shown', false);
        value_ATSAVI.style().set('shown', false);
        BoxC_EVI.style().set('shown', false);
        value_EVI.style().set('shown', false);
        BoxC_EVI2.style().set('shown', false);
        value_EVI2.style().set('shown', false);
        BoxC_OSAVI.style().set('shown', false);
        value_OSAVI.style().set('shown', false);
        BoxC_SAVI.style().set('shown', false);
        value_SAVI.style().set('shown', false);
        BoxC_TSAVI.style().set('shown', false);
        value_TSAVI.style().set('shown', false);
        BoxC_WDRVI.style().set('shown', false);
        value_WDRVI.style().set('shown', false);
    }
  },  
});

//=========================================================***********************************************************************
// CONFIGURACION DE LOS VALORES DE LOS INDICES DE VEGETACIÓN
//==========================================================
//1. EVI (C1, C2 Y L)
var BoxC_ARVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: ARVI', style: St.styleConf }),]
});
var Eps_text= ui.Label({value: 'γ:', style: St.styleTexTConf});
var Eps_in= ui.Textbox({placeholder: '1.0', value: '1.0', style: {maxWidth: '70px'}});
var value_ARVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [Eps_text,Eps_in]
});
//1. EVI (C1, C2 Y L)
var BoxC_EVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: EVI', style: St.styleConf }),]
});
var C1_text= ui.Label({value: 'C1:', style: St.styleTexTConf});
var C1_in= ui.Textbox({placeholder: '6.0', value: '6.0', style: {maxWidth: '70px'}});
var C2_text= ui.Label({value: 'C2:',style: St.styleTexTConf});
var C2_in= ui.Textbox({placeholder: '7.5', value: '7.5', style: {maxWidth: '70px'}});
var L_text= ui.Label({value: 'L:',style: St.styleTexTConf}); 
var L_in= ui.Textbox({placeholder: '1.0', value: '1.0', style: {maxWidth: '70px'}});     
var value_EVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [C1_text,C1_in,C2_text, C2_in,L_text,L_in]
});
//2. EVI2 (C1)
var BoxC_EVI2= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: EVI2', style: St.styleConf }),]
});
var C12_text= ui.Label({value: 'C1:', style: St.styleTexTConf});
var C12_in= ui.Textbox({placeholder: '2.4', value: '2.4', style: {maxWidth: '70px'}});
var value_EVI2= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [C12_text,C12_in]
});
//3. OSAVI (X)
var BoxC_OSAVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: OSAVI', style: St.styleConf }),]
});
var X_text= ui.Label({value: 'X:', style: St.styleTexTConf});
var X_in= ui.Textbox({placeholder: '0.16', value: '0.16', style: {maxWidth: '70px'}});
var value_OSAVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [X_text,X_in]
});
//4. SAVI (L)
var BoxC_SAVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: SAVI', style: St.styleConf }),]
});
var LS_text= ui.Label({value: 'L:', style: St.styleTexTConf});
var LS_in= ui.Textbox({placeholder: '0.5', value: '0.5', style: {maxWidth: '70px'}});
var value_SAVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [LS_text,LS_in]
});
//5.ATSAVI (a, b, X)
var BoxC_ATSAVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: ATSAVI', style: St.styleConf }),]
});
var a_text= ui.Label({value: 'a:', style: St.styleTexTConf});
var a_in= ui.Textbox({placeholder: '1', value: '1', style: {maxWidth: '70px'}});
var b_text= ui.Label({value: 'b:',style: St.styleTexTConf});
var b_in= ui.Textbox({placeholder: '0', value: '0', style: {maxWidth: '70px'}});
var XS_text= ui.Label({value: 'X:',style: St.styleTexTConf}); 
var XS_in= ui.Textbox({placeholder: '0.08', value: '0.08', style: {maxWidth: '70px'}});     
var value_ATSAVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [a_text,a_in,b_text, b_in,XS_text,XS_in]
});
//5.TSAVI (a, b)
var BoxC_TSAVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: TSAVI', style: St.styleConf }),]
});
var aT_text= ui.Label({value: 'a:', style: St.styleTexTConf});
var aT_in= ui.Textbox({placeholder: '1.0', value: '1.0', style: {maxWidth: '70px'}});
var bT_text= ui.Label({value: 'b:',style: St.styleTexTConf});
var bT_in= ui.Textbox({placeholder: '0.0', value: '0.0', style: {maxWidth: '70px'}});
var value_TSAVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [aT_text,aT_in,bT_text, bT_in]
});
//6. WDRVI (alpha)
var BoxC_WDRVI= ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {backgroundColor: "#2c8aaf", color: "#ffffff", margin: '2px 5px'},//border: '0.5px solid #000000' 
  widgets:
    [ui.Label({ value:'Settings: WDRVI', style: St.styleConf }),]
});
var alpha_text= ui.Label({value: '∝:', style: St.styleTexTConf});
var alpha_in= ui.Textbox({placeholder: '0.2', value: '0.2', style: {maxWidth: '70px'}});
var value_WDRVI= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [alpha_text,alpha_in]
});
// archivo vector desde GEE*******************************************************************************

var VectorFileUse=ui.Checkbox({
  label:'Use vector file from GEE',
  value: false,
  style: {color: '#dc7a55', width: '200px'},
  onChange: function(value) {
    if (value === true) { //es para el ARVI
        value_VF.style().set('shown', true);
    }
    else if  (value === false) { //es para el ATSAVI
        value_VF.style().set('shown', false);
    }
  }
})
var VF_text= ui.Label({value: 'URL GEE:', style: St.styleTexTConf});
var VF_in= ui.Textbox({placeholder: 'projects/calcium-verbena-328905/assets/Bate', value: 'projects/calcium-verbena-328905/assets/Bate', style: {maxWidth: '200px'}});
var value_VF= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [VF_text,VF_in]
});
//7. LinearRegresion (alpha)******************************************************************************************************
var Regresion=ui.Checkbox({label:'regression map',value: false,style: {color: '#dc7a55', width: '110px'}})
var RegressionSelect = ui.Select({ 
  items: [ 
          { label:"1. Lineal (Y = a*VI + b)", value: 0 }, 
          { label:"2. Quadratic (Y= a + b*VI + c*VI^2)", value: 1}, 
          { label:"3. Potential (Y= a * VI^b)", value: 2}, 
          { label:"4. Exponential (Y= a*Exp^(b*VI)", value: 3}, 
  ],
  placeholder : 'Selecciona un modelo de regresion',
  value: 0,
  style: {width: '200px', border: '1px solid darkgray'},

}); 
var RegresionT= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {position: 'top-center'},
      widgets: [Regresion, RegressionSelect]
  });
 

var Ra_text= ui.Label({value: 'a:', style: St.styleTexTConf});
var Ra_in= ui.Textbox({placeholder: '1.0', value: '1.0', style: {maxWidth: '70px'}});
var Rb_text= ui.Label({value: 'b:',style: St.styleTexTConf});
var Rb_in= ui.Textbox({placeholder: '0.0', value: '0.0', style: {maxWidth: '70px'}});
var value_Regresion= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [Ra_text,Ra_in,Rb_text, Rb_in]
}); 
var Rac_text= ui.Label({value: 'a:', style: St.styleTexTConf});
var Rac_in= ui.Textbox({placeholder: '0.0', value: '0.0', style: {maxWidth: '70px'}});
var Rbc_text= ui.Label({value: 'b:',style: St.styleTexTConf});
var Rbc_in= ui.Textbox({placeholder: '1.0', value: '1.0', style: {maxWidth: '70px'}});
var Rcc_text= ui.Label({value: 'c:',style: St.styleTexTConf});
var Rcc_in= ui.Textbox({placeholder: '0.0', value: '0.0', style: {maxWidth: '70px'}});
var value_RegresionC= ui.Panel({
   layout: ui.Panel.Layout.flow('horizontal'),
   style: {backgroundColor: "#eceaea", margin: '2px 5px'},
   widgets: [Rac_text,Rac_in,Rbc_text, Rbc_in,Rcc_text, Rcc_in]
}); 

Regresion.onChange(function(value) {
    if (value===true) {
        if (RegressionSelect.getValue()===1){
            value_Regresion.style().set('shown', false);
            value_RegresionC.style().set('shown', true);
        }
        else {
            value_Regresion.style().set('shown', true);
            value_RegresionC.style().set('shown', false); 
        }         
        
        RegressionSelect.onChange(function(value2) {
            if (value2 === 1) { 
                value_Regresion.style().set('shown', false);
                value_RegresionC.style().set('shown', true);        
            }
            else   { 
                value_Regresion.style().set('shown', true);
                value_RegresionC.style().set('shown', false);    
          }  
      })      
    }
    else {
        value_Regresion.style().set('shown', false);
        value_RegresionC.style().set('shown', false); 
        RegressionSelect.onChange(function(value2) {
            if (value2 < 3) { 
                    value_Regresion.style().set('shown', false);
                    value_RegresionC.style().set('shown', false); 
           }
        })
    }
});
//=====================================*******************************************************************************************
//SE CREAN LOS TEXTBOX DE MAX Y MIN
//====================================
var minBox_C1 =ui.Textbox({value:0, placeholder: '0', style:St.TexBoxValueIndex});
var maxBox_C1 =ui.Textbox({value:1, placeholder: '1', style:St.TexBoxValueIndex});
var DownloadC1_Button = ui.Label({value:"⇓", style: St.ButonDonwload});
var DownloadC3_Button = ui.Label({value:"⇓", style:St.ButonDonwload});
var DownloadC2_Button = ui.Label({value:"⇓", style:St.ButonDonwload});
var valC1= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {backgroundColor: "#eceaea", position: 'top-center',shown: false},
      widgets: [ui.Label('VI map display',{backgroundColor: "#eceaea",color: '#000000'}),ui.Label(' min:',{backgroundColor: "#eceaea",color: '#615b58'}),minBox_C1,ui.Label('max:',{backgroundColor: "#eceaea",color: '#615b58'}), maxBox_C1]
  });
var IndiceHelp= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {position: 'top-center'},
      widgets: [indice,VIHelp_Button]
  });
var SensorT= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {position: 'top-center'},
      widgets: [Sensor,SATHelp_Button]
  });
var Pr=ui.Checkbox({label:'Calculate the weighting factor',value: false,style: {color: '#dc7a55', width: '270px'}})
var PrT= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {position: 'top-center'},
      widgets: [Pr,PRHelp_Button]
  });  
var Pol=ui.Checkbox({label:'Filter images covering the entire polygon',value: true,style: {color: '#dc7a55', width: '280px'}})
var PolT= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {position: 'top-center'},
      widgets: [Pol]
  });
var downloadVFile= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {backgroundColor: "#eceaea", position: 'top-center',shown: true},
      widgets: [ui.Label('download vector file',{backgroundColor: "#eceaea", color: '#615b58',width: '270px'}),DownloadC2_Button]
  });
var downloadIV= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {backgroundColor: "#eceaea", position: 'top-center',shown: true},
      widgets: [ui.Label('download vegetation index map',{backgroundColor: "#eceaea", color: '#615b58',width: '270px'}),DownloadC1_Button]
  });
var downloadRM= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {backgroundColor: "#eceaea", position: 'top-center',shown: true},
      widgets: [ui.Label('download Regression map',{backgroundColor: "#eceaea", color: '#615b58',width: '270px'}), DownloadC3_Button]
  });
//seleccion de indice y sensor   
panel.add(indiceT);
panel.add(SensorT);
panel.add(IndiceHelp);
//AJUSTES
panel.add(BoxC_ARVI);
panel.add(value_ARVI);
panel.add(BoxC_ATSAVI);
panel.add(value_ATSAVI);
panel.add(BoxC_EVI);
panel.add(value_EVI);
panel.add(BoxC_EVI2);
panel.add(value_EVI2);
panel.add(BoxC_OSAVI);
panel.add(value_OSAVI);
panel.add(BoxC_SAVI);
panel.add(value_SAVI);
panel.add(BoxC_TSAVI);
panel.add(value_TSAVI);
panel.add(BoxC_WDRVI);
panel.add(value_WDRVI);
panel.add(VectorFileUse);
panel.add(value_VF);
panel.add(RegresionT);
panel.add(value_Regresion);
panel.add(value_RegresionC);
BoxC_ARVI.style().set('shown', false);
value_ARVI.style().set('shown', false);
BoxC_ATSAVI.style().set('shown', false);
value_ATSAVI.style().set('shown', false);
BoxC_EVI.style().set('shown', false);
value_EVI.style().set('shown', false);
BoxC_EVI2.style().set('shown', false);
value_EVI2.style().set('shown', false);
BoxC_OSAVI.style().set('shown', false);
value_OSAVI.style().set('shown', false);
BoxC_SAVI.style().set('shown', false);
value_SAVI.style().set('shown', false);
BoxC_TSAVI.style().set('shown', false);
value_TSAVI.style().set('shown', false);
BoxC_WDRVI.style().set('shown', false);
value_WDRVI.style().set('shown', false);
value_Regresion.style().set('shown', false);
value_RegresionC.style().set('shown', false);
value_VF.style().set('shown', false);
panel.add(Pol);
panel.add(PrT);
panel.add(valC1);
panel.add(downloadIV);
panel.add(downloadVFile);
panel.add(downloadRM);
downloadVFile.style().set('shown', false);
downloadIV.style().set('shown', false);
downloadRM.style().set('shown', false);
// 5.4. Crea botón para procesar, con la función de llamada Procesa y lo agrega al panel principal*****************
var proc = ui.Button('Calculate', Procesa, false, {margin: 'auto', padding: '5px', width: '110px', fontSize: '15px'});
var proc2 = ui.Button('New polygon', Nuevo, false, {margin: 'auto', padding: '5px', width: '110px', fontSize: '15px'});
var proc3 = ui.Button('Edit polygon', Edita, false, {margin: 'auto', padding: '5px', width: '110px', fontSize: '15px'});
proc3.setDisabled(true);
proc2.setDisabled(true);
var procHelp= ui.Panel({
      layout: ui.Panel.Layout.flow("horizontal"),
      style: {position: 'top-center'},
      widgets: [proc2,proc3,proc]
  });
panel.add(procHelp);


//FUNCION SI CAMBIA LA IMAGEN  (items,uris,i,ivVis,nombreInd,recorte,lis_Indice,lis_S2A )****************************************
function buildImageSelect(items,uris,i,ivVis, reVis, nombreInd,recorte,Pr_Response,lis_S2A,prodVis, capa4, indSel,scala,lis_Regresion, Re_Response, nombreRegresion){ 
       var imageSelect = ui.Select({ 
           items: items, 
           value:0,
           placeholder: "Lista de Imagenes", 
           style: { width: '95%', margin: '10px 10px 10px 10px', position: 'top-center'},
           onChange: function (value){
                    var iv=ee.Image(uris.get(value));
                    var ivR=ee.Image(lis_Regresion.get(value));
                    var capa1=ui.Map.Layer(iv.clip(recorte), ivVis, nombreInd);
                    var capa2= ui.Map.Layer((ee.Image(lis_S2A.get(value))).clip(recorte), St.vizParams, 'RGB');
                    var prodrel = imp2.ProdRelativa(iv, recorte,nombreInd);
                    var capa3=ui.Map.Layer(ee.Image(prodrel), prodVis, 'Weighting Factor');  
                    var capa5=ui.Map.Layer(ivR.select(nombreRegresion).clip(recorte), reVis, nombreRegresion);
                    if (Pr_Response===true) {
                        if (Re_Response===true) { 
                              Map.layers().reset([capa2,capa1, capa3,capa5, capa4]); 
                        }
                        else {
                              Map.layers().reset([capa2,capa1, capa3, capa4]);
                        }
                    }
                    else {
                        if (Re_Response===true) { 
                              Map.layers().reset([capa2, capa1,capa5, capa4]);  
                        }
                        else {
                              Map.layers().reset([capa2,capa1, capa4]);
                        }
                    }    
                    DownloadC1_Button.setUrl(getURLNDVI(iv,nombreInd, recorte,scala,nombreInd ));
                    DownloadC3_Button.setUrl(getURLNDVI(ee.Image(lis_Regresion.get(value)),nombreRegresion, recorte,scala,nombreRegresion));
                    minBox_C1.setDisabled(false);
                    maxBox_C1.setDisabled(false); 
              //=============================================
              //CUANDO CAMBIA LA CLASE 1
              //============================================
              minBox_C1.onChange(function(value) {
      
                    ivVis.min =parseFloat(value);
                    capa1.setVisParams(ivVis);
                    var pan1 = Leyenda(nombreInd, ivVis);
                    panL1.widgets().set(2, pan1.widgets().get(2));
                    capa1 = ui.Map.Layer(iv.clip(recorte), ivVis, nombreInd);
                    if (Pr_Response===true) {
                        if (Re_Response===true) { 
                              Map.layers().reset([capa2,capa1, capa3,capa5, capa4]); 
                        }
                        else {
                              Map.layers().reset([capa2,capa1, capa3, capa4]);
                        }
                    }
                    else {
                        if (Re_Response===true) { 
                              Map.layers().reset([capa2, capa1,capa5, capa4]);  
                        }
                        else {
                              Map.layers().reset([capa2,capa1, capa4]);
                        }         
                    }
              });
              maxBox_C1.onChange(function(value) {
                    ivVis.max =parseFloat(value);
                    capa1.setVisParams(ivVis);
                    var pan1 = Leyenda(nombreInd, ivVis);
                    panL1.widgets().set(2, pan1.widgets().get(2));
                    capa1 = ui.Map.Layer(iv.clip(recorte), ivVis, nombreInd);
                    if (Pr_Response===true) {
                        if (Re_Response===true) { 
                              Map.layers().reset([capa2,capa1, capa3,capa5, capa4]); 
                        }
                        else {
                              Map.layers().reset([capa2,capa1, capa3, capa4]);
                        }
                    }
                    else {
                        if (Re_Response===true) { 
                              Map.layers().reset([capa2, capa1,capa5, capa4]);  
                        }
                        else {
                              Map.layers().reset([capa2,capa1, capa4]);
                        }         
                    }
              });  
            }
       });
       return imageSelect; 
}

function buildSelectDisplay(number) {
  if (number==1) {
    minBox_C1.setValue(0);
    maxBox_C1.setValue(1);
  }
  else if (number==2) {
    minBox_C1.setValue(-1);
    maxBox_C1.setValue(0);
  }
  else if (number==3) {
    minBox_C1.setValue(0);
    maxBox_C1.setValue(2);
  }
  else if (number==4) {
    minBox_C1.setValue(-1);
    maxBox_C1.setValue(1);
  }
  else if (number==5) {
    minBox_C1.setValue(0);
    maxBox_C1.setValue(100);
  }
  else if (number==6) {
    minBox_C1.setValue(0);
    maxBox_C1.setValue(40);
  }
  else if (number==7) {
    minBox_C1.setValue(-0.5);
    maxBox_C1.setValue(1);
  }
   valC1.style().set('shown', true); //clase 1
   minBox_C1.setDisabled(true);
   maxBox_C1.setDisabled(true);
}

// ===============================================================================**************************************
// 6. LAMADA DE FUNCIONES
// ===============================================================================************************************

var fecha1;
var fecha2;
var ivs;
var nombreInd;
var panL;
var panL1;
var prodVis;
var ivVis;
var capa0;

//=============================================================********************************************************
function Nuevo()
{
  Map.clear(); 
  //Map.setControlVisibility(true, true, true, true, true, true, true);
  drawingTools.clear();
  proc3.setDisabled(true)
  Map.setCenter(-101.59, 24.10, 6);
  valC1.style().set('shown',false);
  var inifap_leyend = ui.Label('INIFAP');
  inifap_leyend.style().set('position', 'bottom-left');
  inifap_leyend.setUrl('https://www.gob.mx/inifap');
  
  Map.add(inifap_leyend);
  var titulo = ui.Label({value:'Draw the polygons of the agricultural area using the drawing tools', style:{color:'red'}});
titulo.style().set('position', 'bottom-center');
Map.add(titulo);
 Map.style().set('cursor', 'hand');
}

function Edita()//****************************************************************************************************************
{
  Map.clear(); 
  //Map.setControlVisibility(true, true, true, true, true, true, true);
  //drawingTools.clear();
 // Map.setCenter(-101.59, 24.10, 6);
  valC1.style().set('shown',false);
    Map.drawingTools().layers().forEach(function(layer) {
    layer.setShown(true);
  });
  var inifap_leyend = ui.Label('INIFAP');
  inifap_leyend.style().set('position', 'bottom-left');
  inifap_leyend.setUrl('https://www.gob.mx/inifap');
  
  Map.add(inifap_leyend);
  var titulo = ui.Label({value:'Draw the polygons of the agricultural area using the drawing tools', style:{color:'red'}});
titulo.style().set('position', 'bottom-center');
Map.add(titulo);
 Map.style().set('cursor', 'hand');
}

function Procesa()//*********************************************************************************************************
{
  var panel1 = ui.Panel([], ui.Panel.Layout.Flow('vertical'),{position:'top-right'});
  
  proc2.setDisabled(false)
  // 1) parar el dibujo e importar la capa creada por el usuario
  drawingTools.stop();
  var modulo = drawingTools.getMap();
  var om= drawingTools.toFeatureCollection(modulo);
  //Ocultar la capa creada por el usuario
  Map.drawingTools().layers().forEach(function(layer) {
    layer.setShown(false);
  });
  
  //comenzamos a importar todo lo relacionado a las imagenes
  var umbral = parseFloat(kh_in.getValue());
  var fecha = [];
  fecha1 = fechain.getValue();
  fecha2 = fechafin.getValue();
  fecha.push(fecha1); 
  fecha.push(fecha2);
  var indSel = indice.getValue();
  var SenSel = Sensor.getValue();
  var Pr_Response=Pr.getValue();
  var Re_Response=Regresion.getValue();
  var ReModel_Response=RegressionSelect.getValue();
  var Pol_Response=Pol.getValue();
  var VF_Response=VectorFileUse.getValue();
  var recorte;
  if (VF_Response===true){
    recorte=ee.FeatureCollection (VF_in.getValue());
    proc3.setDisabled(true)
  }
  else {
     recorte = (om);
     proc3.setDisabled(false)
  }
  
  var S2B;
  var S2A;
  var scala;
  var VistaRGB;
  // Crea colección Sentinel 2 filtrado por fecha, área de interés y umbral de nubes=============
  if (SenSel===0) {
     S2B = ColeccionImagenSR(fecha, recorte, umbral);
     scala=30;
  }
  else if  (SenSel==1) {
     S2B = ColeccionImagenSentinelSR (fecha, recorte, umbral);
     scala=10;
  }
  else if (SenSel==2) {
     S2B = ColeccionImagenAMBOS(fecha, recorte, umbral);
     scala=10;     
  }  
  else {
     S2B = ColeccionImagenSR_4(fecha, recorte, umbral);
     scala=30;
  }
  print(S2B)
  //filtra imagenes que cubren todo el poligono
  if  (Pol_Response===false) {
        S2B=S2B;
  }
  else {
        S2B=ee.ImageCollection(S2B.filter(ee.Filter.contains('.geo', recorte.geometry())));
  }
  var lis_S2A=S2B.toList(S2B.size());
  S2A=ee.ImageCollection(lis_S2A)
  var im = S2A.first();
  var X_val;
  var alpha_val;
  var a_val;
  var b_val;
  //Calculamos los indices de vegetación
  if (indSel === 0) { //Atmospherically ResistanA GLOBAL CALCULATOR OF VEGETATION INDICES t Vegetation Index-ARVI
    ivs = ee.ImageCollection(S2A.map(ARVI));
    nombreInd = 'ARVI';
    buildSelectDisplay(7);
  }
  else if (indSel == 1) { //Enhanced Vegetation index-EVI
    ivs = ee.ImageCollection(S2A.map(ATSAVI));
    nombreInd = 'ATSAVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 2) { //Enhanced Vegetation index-EVI
    ivs = ee.ImageCollection(S2A.map(imp2.DVI));
    nombreInd = 'DVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 3) { //Enhanced Vegetation index-DVI
    ivs = ee.ImageCollection(S2A.map(EVI));
    nombreInd = 'EVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 4) { //Enhanced Vegetation index-EVI2
    ivs = ee.ImageCollection(S2A.map(EVI2));
    nombreInd = 'EVI2';
    buildSelectDisplay(1);
  }
  else if (indSel == 5) { //Green Normalized Difference Vegetation Index-GNDVI
    ivs = ee.ImageCollection(S2A.map(imp2.GNDVI));
    nombreInd = 'GNDVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 6) { //Modified Soil Adjusted Vegetation Index-MSAVI2
    ivs = ee.ImageCollection(S2A.map(imp2.MSAVI2));
    nombreInd = 'MSAVI2';
    buildSelectDisplay(1);  
  }
  else if (indSel == 7) { //Moisture Stress Index-MSI
    ivs = ee.ImageCollection(S2A.map(imp2.MSI));
    nombreInd = 'MSI';
    buildSelectDisplay(3); 
  }
  else if (indSel == 8) { // Modiﬁed triangular vegetation index-MTVI
    ivs = ee.ImageCollection(S2A.map(imp2.MTVI));
    nombreInd = 'MTVI';
    buildSelectDisplay(7); 
  }
  else if (indSel == 9) { //Modiﬁed Triangular Vegetation Index-MTVI2
    ivs = ee.ImageCollection(S2A.map(imp2.MTVI2));
    nombreInd = 'MTVI2';
    buildSelectDisplay(7);   
  }
  else if (indSel == 10) { //Normalized Difference Tillage Index (NDTI)
    ivs = ee.ImageCollection(S2A.map(imp2.NDTI));
    nombreInd = 'NDTI';
    buildSelectDisplay(1);
  }
  else if (indSel == 11) { //Normalized Difference Vegetation Index- NDVI
    ivs = ee.ImageCollection(S2A.map(imp2.NDVI));
    nombreInd = 'NDVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 12) { //Normalized Difference Water Index-NDWI
    ivs = ee.ImageCollection(S2A.map(imp2.NDWI));
    nombreInd = 'NDWI';
    buildSelectDisplay(7);
  }
  else if (indSel == 13) { //Optimized Soil Adjusted Vegetation Index-OSAVI
    ivs = ee.ImageCollection(S2A.map(OSAVI));
    nombreInd = 'OSAVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 14) { //Optimized Soil Adjusted Vegetation Index-OSAVI
    ivs = ee.ImageCollection(S2A.map(imp2.RDVI));
    nombreInd = 'RDVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 15) { //Redness Index-RI
    ivs = ee.ImageCollection(S2A.map(imp2.RI));
    nombreInd = 'RI';
    buildSelectDisplay(4);
  }
  else if (indSel == 16) { //Ratio Vegetation Index-RVI
    ivs = ee.ImageCollection(S2A.map(imp2.RVI));
    nombreInd = 'RVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 17) { //Soil Adjusted Vegetation Index-SAVI
    ivs = ee.ImageCollection(S2A.map(SAVI));
    nombreInd = 'SAVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 18) { //Soil Adjusted Vegetation Index-SAVI
    ivs = ee.ImageCollection(S2A.map(imp2.TVI));
    nombreInd = 'TVI';
    buildSelectDisplay(6);
  }
  else if (indSel == 19) { //Transformed Soil Adjusted Vegetation Index-TSAVI
    ivs = ee.ImageCollection(S2A.map(TSAVI));
    nombreInd = 'TSAVI';
    buildSelectDisplay(1);
  }
  else if (indSel == 20) { //Visible Atmospherically Resistant Index-VARI
    ivs = ee.ImageCollection(S2A.map(imp2.VARI));
    nombreInd = 'VARI';
    buildSelectDisplay(4);
  }
  else if (indSel == 22) { //Vegetation Index Number-VIN
    ivs = ee.ImageCollection(S2A.map(imp2.VIN));
    nombreInd = 'VIN';
    buildSelectDisplay(6);
  }
  else if (indSel == 23) { //Wide Dynamic Range Vegetation Index-WDRVI
    ivs = ee.ImageCollection(S2A.map(WDRVI));
    nombreInd = 'WDRVI';
    buildSelectDisplay(4);
  }
  
  //se agrega el link an boton de descarga y se recuperan valores de max y min
  downloadIV.style().set('shown', true);
  downloadVFile.style().set('shown', true);  
  minPalete=parseFloat(minBox_C1.getValue());
  maxPalete=parseFloat(maxBox_C1.getValue());
    
  // Dibuja parcelas como imagen
  var parcelas = ee.Image().byte();
  bordes = parcelas.paint({
    featureCollection: recorte,
    color: 1,
    width: 1
  });
  Map.clear(); 
  Map.centerObject(recorte, 13);
  Map.setControlVisibility(true, true, false, true, true, true, false);
  // Crea título del Mapa.
  var titulo = ui.Label('Click inside a polygon to plot the time series');
  titulo.style().set('position', 'top-center');
  Map.add(titulo);
  Map.style().set('cursor', 'crosshair');
  var mosf = ee.Date(ee.Image(S2A.first()).get('system:time_start'));
  var iv = ivs.first();
  im = S2A.filterDate(mosf.advance(-0.5, 'hour'), mosf.advance(0.5, 'hour'));
  var prodrel = imp2.ProdRelativa(iv, recorte,nombreInd); //Productividad
  //se obtiene el mapa de modelo de regresion
  var RegresionValue
  var nombreRegresion
  var MinReg
  var MaxReg
  var Rega=parseFloat(Ra_in.getValue())
  var Regb=parseFloat(Rb_in.getValue())
  var Regac=parseFloat(Rac_in.getValue())
  var Regbc=parseFloat(Rbc_in.getValue())
  var Regcc=parseFloat(Rcc_in.getValue())
  var minRV=ee.Number(minPalete)
  var maxRV=ee.Number(maxPalete)
  
  if (ReModel_Response===0){ 
      RegresionValue=ee.ImageCollection(ivs.map(RegresionF)); //regresion 
      nombreRegresion='linealModel';
      MinReg=minRV.multiply(Rega).add(Regb);
      MaxReg=maxRV.multiply(Rega).add(Regb);
  }
  else if (ReModel_Response==1){ 
      RegresionValue=ee.ImageCollection(ivs.map(RegresionFc)); //regresion 
      nombreRegresion='QuadraticModel';
      MinReg=minRV.multiply(Regbc).add(Regac).add((minRV.pow(2)).multiply(Regcc));
      MaxReg=maxRV.multiply(Regbc).add(Regac).add((maxRV.pow(2)).multiply(Regcc));
  }
  else if (ReModel_Response==2){ 
      RegresionValue=ee.ImageCollection(ivs.map(RegresionFp)); //regresion 
      nombreRegresion='PotentialModel';  
      if (minRV<=ee.Number(0)) {
        MinReg=ee.Number(0).multiply(minRV);
      }
      else {
        MinReg=(minRV.pow(Regb)).multiply(Rega);
      }
      MaxReg=(maxRV.pow(Regb)).multiply(Rega);
  }
  else if (ReModel_Response==3){ 
      RegresionValue=ee.ImageCollection(ivs.map(RegresionFe)); //regresion 
      nombreRegresion='ExponentialModel';  
      MinReg=((minRV.multiply(Regb)).exp()).multiply(Rega)
      MaxReg=((maxRV.multiply(Regb)).exp()).multiply(Rega)
  } 
  
  var Regresionfirst=RegresionValue.first();

  //Se agregan las cuatro capas al mapa 
  var tituloPR = 'Weighting Factor (adim.)';
  var tituloPR2 = 'regression values';
  var Capa5=3;
  var reVis;
  var panL2;
  prodVis = {min : 0.5, max : 1.5, palette : St.paletaProd};
  ivVis = {min :minPalete, max : maxPalete, palette : St.paletaIV};
  Map.addLayer(im.mosaic().clip(recorte), St.vizParams, 'RGB'); //RGB
  Map.addLayer(iv.clip(recorte), ivVis, nombreInd); //Indice


  if (Pr_Response===true) {
    Map.addLayer(ee.Image(prodrel), prodVis, 'Weighting Factor'); //Productividad Relativa
    panL = Leyenda(tituloPR, prodVis);
    Capa5=4;
    Map.add(panL);
    if (Re_Response===true) { 
        reVis = {min :MinReg.getInfo(), max : MaxReg.getInfo(), palette : St.ETO};  
        Map.addLayer(Regresionfirst.select(nombreRegresion).clip(recorte), reVis, nombreRegresion); //Productividad Relativa
        panL2 = Leyenda(tituloPR2, reVis);
        Capa5=5;
        Map.add(panL2); 
        downloadRM.style().set('shown', true);

    }
    else {
        downloadRM.style().set('shown', false);
    }
  }
  else {
    if (Re_Response===true) { 
        reVis = {min :MinReg.getInfo(), max : MaxReg.getInfo(), palette :St.ETO};  
        Map.addLayer(Regresionfirst.select(nombreRegresion).clip(recorte), reVis, nombreRegresion); //Productividad Relativa
        panL2 = Leyenda(tituloPR2, reVis);
        Capa5=4;
        Map.add(panL2); 
        downloadRM.style().set('shown', true);

    }
    else {
        downloadRM.style().set('shown', false);
    }
  }
  
  Map.addLayer(bordes, {palette: 'black'}, 'Polygons');
  panL1 = Leyenda(nombreInd, ivVis);
  panL = Leyenda(tituloPR, prodVis);
  Map.add(panL1);
  var capa4=ui.Map.Layer(bordes, {palette: 'black'}, 'Polygon'); //poligono
  //se generan las listas para extraer informacion
  var uris=ivs.toList(ivs.size());
  var lis_Regresion=RegresionValue.toList(RegresionValue.size());
  var imageSelect = null;
  S2A.evaluate(function(col){ //se agregan las fechas a un listbox
            var items = []; 
            var i=0;
            col.features.forEach(function(feature){ 
                  var label =i +". " +feature.properties.fecha + " / " 
                   +feature.properties.sensor + " / " + feature.properties.nubosidad + " % ";
                  var value = i;
               items.push({label: label, value:value});
              i=i+1;
            }); 
            imageSelect = buildImageSelect(items,uris,i,ivVis,reVis,nombreInd,recorte,Pr_Response,lis_S2A,prodVis, capa4,indSel,scala, lis_Regresion,Re_Response, nombreRegresion);
            var image2Label = ui.Label('Select a date to view VI map');
            var query2Panel = ui.Panel();
            panel1.add(image2Label)
            panel1.add(imageSelect)
            panel1.add(query2Panel)
  });  
  DownloadC1_Button.setUrl(getURLNDVI(iv,nombreInd, recorte,scala));
  DownloadC3_Button.setUrl(getURLNDVI(Regresionfirst,nombreRegresion, recorte,scala));
  DownloadC2_Button.setUrl(getVectorFile(recorte));
  Map.onClick(function(coords) 
  {
    pangr.style().set('shown', true);
    punto = ee.Geometry.Point(coords.lon, coords.lat);
    //var sergioomar=new ee.FeatureCollection([geometry])
    var filtered = recorte.filterBounds(punto);
    var area2=filtered.geometry().area().divide(10000 );
    var perimetro=filtered.geometry().perimeter().divide(1000 );
    var Eleva = elevation.reduceRegion({
              reducer: ee.Reducer.mean(),
              geometry: punto,
              scale: 100
    });
    var ElevaValue = Eleva.get('elevation');
    
    Map.layers().set(Capa5, filtered.style(St.HIGHLIGHT_STYLE)); // SE AGREGA CAPA DE SELECCION 
    var GrafIV;
        
    var reducers=  ee.Reducer.mean().combine({
        reducer2: ee.Reducer.stdDev(),
        sharedInputs: true
    });
    
    if (Re_Response===true) { 
          // Crea gráfica de índices de vegetación
        GrafIV = ui.Chart.image.series(RegresionValue.select([nombreInd, nombreRegresion]),filtered,ee.Reducer.mean(), 10)
            .setOptions({
              title: 'vegetation index time series',
              series:{
                    0: {targetAxisIndex: 0, curveType: 'function', pointSize: 4,lineWidth: 1, }, //Etapa (Derecho)
                    1: {targetAxisIndex: 1, curveType: 'function', pointSize: 4, lineWidth: 1}  // NDVI (Izquierdo)
                     },              
              hAxis: {title: 'Date', format: 'MM-yyyy'},
              vAxes: {
                0: {
                  title: nombreInd,
                  baseline: 0,
                  titleTextStyle: {italic: false, bold: true, color: '11591a'},
                  viewWindow: {max: maxPalete}
                },
                1: {
                  title: 'Regression value',
                  baseline: 0,
                  titleTextStyle: {italic: false, bold: true, color: '115075'},
                  viewWindow: {max: MaxReg.getInfo()}
                },
              },
              chartArea: {backgroundColor: 'EBEBEB'},
              interpolateNulls: true,
          });
    }
    else {
           // Crea gráfica de índices de vegetación
        GrafIV = ui.Chart.image.series(ivs,filtered,reducers, 10)
            .setOptions({
              title: 'vegetation index time series',
              hAxis: {title: 'Date', format: 'MM-yyyy'},
              vAxis: {title: nombreInd, baseline: 0, viewWindow: { max: maxPalete}}, //
              pointSize: 4,
              curveType: 'function',
              lineWidth: 1,
              chartArea: {backgroundColor: 'EBEBEB'},
              interpolateNulls: true,
          });
    }
    pangr.widgets().set(0,(ui.Button({
                label: 'Close',
                onClick: function() {
                  pangr.style().set('shown', false);
                }
              })));    
    //pangr.widgets().set(1, ui.Label({value:'Click on a point on the graph to see the image of that date', style:{color:'green'}}));
    pangr.widgets().set(1, GrafIV);
    ElevaValue.evaluate(function(result2) {
    area2.evaluate(function(result){
    pangr.widgets().set(2, ui.Label({ 
                  value: 'Area(Ha): ' + result.toFixed(2)+';  '+ '  Elevation(msl): ' + result2.toFixed(2),
                  style: {stretch: 'vertical'}
                }));
      })});
      
  });
  Map.add(panel1);
  Map.add(pangr);
}
var inifap_leyend = ui.Label('INIFAP');
inifap_leyend.style().set('position', 'bottom-left');
inifap_leyend.setUrl('https://www.gob.mx/inifap')
var titulo = ui.Label({value:'Draw the polygons of the agricultural area using the drawing tools', style:{color:'red'}});
titulo.style().set('position', 'bottom-center');
Map.add(titulo);
Map.add(inifap_leyend);
Map.setCenter(-101.59, 24.10, 6);
drawingTools.setDrawModes(['polygon']);

ui.root.insert(0, panel);
