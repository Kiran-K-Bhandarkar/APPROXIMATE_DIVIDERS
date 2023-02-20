LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY INEXACT_CELL_Q IS PORT
(
	X, Y, BIN: IN STD_LOGIC;
	D, BOUT	 : OUT STD_LOGIC
);
END ENTITY INEXACT_CELL_Q;

ARCHITECTURE BEHAVIORAL OF INEXACT_CELL_Q IS
BEGIN
	BOUT <= NOT(X);
	D    <= '0';
END ARCHITECTURE BEHAVIORAL;

