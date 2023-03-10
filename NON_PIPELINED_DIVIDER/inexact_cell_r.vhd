LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY INEXACT_CELL_R IS PORT
(
	QS, X, D : IN STD_LOGIC;
	R 	     : OUT STD_LOGIC
);
END ENTITY INEXACT_CELL_R;

ARCHITECTURE BEHAVIORAL OF INEXACT_CELL_R IS
BEGIN
	R <= QS XOR X;
END ARCHITECTURE BEHAVIORAL;


