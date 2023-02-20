LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY EXACT_CELL_R IS PORT
(
	QS, X, D : IN STD_LOGIC;
	R 	     : OUT STD_LOGIC
);
END ENTITY EXACT_CELL_R;

ARCHITECTURE BEHAVIORAL OF EXACT_CELL_R IS
BEGIN
  OUT_PROCESS: PROCESS(QS, X, D)
   BEGIN
	CASE QS IS 
		WHEN '0' 	=>
		  R <= X;
		WHEN '1' 	=>
		  R <= D;
		WHEN OTHERS =>
		  R <= 'X';
    END CASE;
  END PROCESS OUT_PROCESS;
END ARCHITECTURE BEHAVIORAL;


