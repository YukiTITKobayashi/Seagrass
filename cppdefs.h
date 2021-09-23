/* mod_input.F90 */
/* ----- READ Tide data (m) --------------------------------------------- */
/* #define TIDE_DATA_FILE */
/* #define TIDE_CSV_FILE */
/* #define TIDE_FORMULA1 */
#define TIDE_FORMULA2

/* ----- READ ambient water column & pore water salinity data (psu) ----------------- */
/* #define SAL_CSV_FILE */
#define SAL_FORMULA_FROM_TIDE   /* This preprocessor must be used with TIDE_FORMULA1 or TIDE_FORMULA2 */

/* #define SGFLUX_HOFLER_THODAY */
#define SGFLUX_NATURAL
/* #define SGFLUX_SIMPLE */
