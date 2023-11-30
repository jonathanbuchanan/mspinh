#!/usr/bin/env sage
import argparse
import csv

# Parse arguments.
parser = argparse.ArgumentParser(
  description="Calculate homotopy groups of MSpin^x and write the results to a CSV file. The nth row (starting at zero) is the nth homotopy group. The first value is the number of Z summands, and the second is the number of Z/2Z summands.")
parser.add_argument(
  "--out",
  nargs=1,
  default="pi_MSpinx.csv",
  help="The CSV file to write the resulting data to")
parser.add_argument(
  "n",
  type=int,
  nargs=1,
  help="The number of homotopy groups to compute.")
parser.add_argument(
  "x",
  nargs=1,
  choices=["r", "c", "h"],
  help="The variant of MSpin to compute the homotopy groups of")
args = parser.parse_args()
output_file = args.out[0]
n = args.n[0]
variant = args.x[0]

# Note that we need one extra degree of precision for computing the Poincare
# series of the elephant module.
R.<t> = PowerSeriesRing(ZZ, default_prec=n)

def one_m_tk_inv(k):
  return (1 - t^k).inverse()


# --- Poincare series for Steenrod modules ---
# The Poincare series for A, the mod 2 Steenrod algebra.
def poincare_series_A():
  # The series is the product of (1 - t^{2^k - 1})^{-1} for k >= 1. See Theorem
  # 1.11.b in ABP65.
  result = 1
  k = 1
  while n >= (2^k) - 1:
    result = result * one_m_tk_inv((2^k) - 1)
    k += 1
  return result

# The Poincare series for A/(Sq1, Sq2).
def poincare_series_A_mod_Sq1_Sq2():
  # The series is the product of (1 - t^{2^k - 1})^{-1} * (1 - t^4)^{-1}
  # * (1 - t^6)^{-1} for k >= 3. See Theorem 1.11.c in ABP65.
  result = one_m_tk_inv(4)
  result = result * one_m_tk_inv(6)
  k = 3
  while n >= (2^k) - 1:
    result = result * one_m_tk_inv((2^k) - 1)
    k += 1
  return result

# The Poincare series for A/A(Sq3).
def poincare_series_A_mod_Sq3():
  # The series is the product of (1 - t^{2^k - 1})^{-1} * (1 - t^4)^{-1}
  # * (1 - t^6)^{-1} * f for k >= 3, where f = 1 + t + t^2 + t^3 + t^4. See
  # Theorem 1.11.d in ABP65.
  A_mod_sq1_sq2 = poincare_series_A_mod_Sq1_Sq2()
  f = 1 + t + t^2 + t^3 + t^4
  return A_mod_sq1_sq2 * f

# The Poincare series for A/A(Sq1, Sq3).
def poincare_series_A_mod_Sq1_Sq3():
  # The series is the product of (1 - t^{2^k - 1})^{-1} * (1 - t^2)^{-1}
  # * (1 - t^6)^{-1} for k >= 3.
  # See Stong.
  result = one_m_tk_inv(2)
  result = result * one_m_tk_inv(6)
  k = 3
  while n >= (2^k) - 1:
    result = result * one_m_tk_inv((2^k) - 1)
    k += 1
  return result

# The Poincare series for Q_A aka the upside down question mark aka
# A/A(Sq1, Sq5).
def poincare_series_Q():
  # The series is (1 + t^2 + t^3) * (1 - t^4)^{-1} * (1 - t^6)^{-1} times
  # the product of (1 - t^{2^k - 1})^{-1} for k >= 3.
  f = 1 + t^2 + t^3
  f = f * one_m_tk_inv(4)
  f = f * one_m_tk_inv(6)
  k = 3
  while n >= (2^k) - 1:
    f = f * one_m_tk_inv((2^k) - 1)
    k += 1
  return f

# The Poincare series for E_A aka the elephant aka ASq1 + ASq2.
def poincare_series_E():
  # The series is (1 - t^4)^{-1} * (1 - t^6)^{-1}
  # * (1 + t + 2t^2 + t^3 + t^4 + t^5) times the product of
  # (1 - t^{2^k - 1})^{-1} for k >= 3.
  f = one_m_tk_inv(4)
  f = f * one_m_tk_inv(6)
  f = f * (1 + t + (2 * t^2) + t^3 + t^4 + t^5)
  k = 3
  while n >= (2^k) - 1:
    f = f * one_m_tk_inv((2^k) - 1)
    k += 1
  return f

# The Poincare series for the cohomology of MSpin.
def poincare_series_H_mspin():
  # The series is the product of all (1 - t^k)^{-1} for k > 3 and k - 1 not a
  # power of two. See Theorem 1.11.a in ABP66.
  result = one_m_tk_inv(4)
  for k in range(6, n):
    # Check that k - 1 is not a power of 2
    if ((k - 1) & (k - 2) != 0):
      result = result * one_m_tk_inv(k)
  return result

# The Poincare series for the cohomology of MSpin^c.
def poincare_series_H_mspinc():
  # The series is the product of all (1 - t^k)^{-1} for k > 1 and k - 1 not a
  # power of two.
  result = one_m_tk_inv(2)
  for k in range(4, n):
    # Check that k - 1 is not a power of 2
    if ((k - 1) & (k - 2) != 0):
      result = result * one_m_tk_inv(k)
  return result

# The Poincare series for the cohomology of MSpin^h.
def poincare_series_H_mspinh():
  # The series is the product of all (1 - t^k)^{-1} for k > 1 and k - 1 not a
  # power of two or k = 3.
  result = one_m_tk_inv(2) * one_m_tk_inv(3)
  for k in range(4, n):
    # Check that k - 1 is not a power of 2
    if ((k - 1) & (k - 2) != 0):
      result = result * one_m_tk_inv(k)
  return result


# --- Partitions ---
# The number of partitions of n where 1 does not occurs is equal to the number
# of partitions of n minus the number of partitions of n - 1.
def partitions_without_1(n):
  if n == 0:
    return 1
  else:
    return Partitions(n).cardinality() - Partitions(n - 1).cardinality()


# --- Poincare series for summands splitting MSpin^x ---
# The Poincare series for ko summands in MSpin.
def poincare_series_MSpin_ko():
  # One ko summand, shifted to degree 8k, for each partition of 2k not
  # containing 1.
  result = 0
  for i in range(0, (n // 8) + 1):
    result += partitions_without_1(2 * i) * t^(8 * i)
  return result

# The Poincare series for ko<2> summands in MSpin.
def poincare_series_MSpin_ko2():
  # One ko summand, shifted to degree 8k + 2, for each partition of 2k + 1 not
  # containing 1.
  result = 0
  for i in range(0, (n // 8) + 1):
    result += partitions_without_1((2 * i) + 1) * t^((8 * i) + 2)
  return result

# The Poincare series for HZ/2Z summands in MSpin.
def poincare_series_MSpin_HZ2(ko, ko2):
  # Take the Poincare series for the cohomology of MSpin. Subtract off the
  # cohomology of the ko summands and then divide by the Poincare series of A.
  ZA = (poincare_series_H_mspin()
    - (ko * poincare_series_A_mod_Sq1_Sq2())
    - (ko2 * poincare_series_A_mod_Sq3()))
  Z = ZA * poincare_series_A().inverse()
  return Z

# The Poincare series for ku summands in MSpin^c.
def poincare_series_MSpinc_ku():
  # One ku summand, shifted to degree 4k, for each partition of k.
  result = 0
  for i in range(0, (n // 4) + 1):
    result += Partitions(i).cardinality() * t^(4 * i)
  return result

# The Poincare series for HZ/2Z summands in MSpin^c.
def poincare_series_MSpinc_HZ2(ku):
  # Take the Poincare series for the cohomology of MSpin^c. Subtract off the
  # cohomology of the ku summands and then divide by the Poincare series of A.
  ZA = (poincare_series_H_mspinc()
    - (ku * poincare_series_A_mod_Sq1_Sq3()))
  Z = ZA * poincare_series_A().inverse()
  return Z

# The Poincare series for ksp summands in MSpin^h.
def poincare_series_MSpinh_ksp():
  # One ksp summand, shifted to degree 8k, for each partition of 2k.
  result = 0
  for i in range(0, (n // 8) + 1):
    result += Partitions(2 * i).cardinality() * t^(8 * i)
  return result

# The Poincare series for F summands in MSpin^h.
def poincare_series_MSpinh_F():
  # One F summand, shifted to degree 8k + 4, for each partition of 2k + 1.
  result = 0
  for i in range(0, (n // 8) + 1):
    result += Partitions((2 * i) + 1).cardinality() * t^((8 * i) + 4)
  return result

# The Poincare series for HZ/2Z summands in MSpin^h.
def poincare_series_MSpinh_HZ2(ksp, f):
  # Take the Poincare series for the cohomology of MSpin^h. Subtract off the
  # cohomology of the ksp and F summands and then divide by the Poincare series
  # of A.
  ZA = (poincare_series_H_mspinh()
    - (ksp * poincare_series_Q())
    - (f * poincare_series_E()))
  Z = ZA * poincare_series_A().inverse()
  return Z


# --- Poincare series for homotopy groups ---
# The Poincare series for Z summands in the homotopy groups of ko.
def poincare_series_pi_ko_Z():
  return one_m_tk_inv(4)

# The Poincare series for Z/2Z summands in the homotopy groups of ko.
def poincare_series_pi_ko_Z2():
  return (t * one_m_tk_inv(8)) + (t^2 * one_m_tk_inv(8))

# The Poincare series for Z summands in the homotopy groups of ko<2> (looped
# twice).
def poincare_series_pi_ko2_Z():
  return t^2 * one_m_tk_inv(4)

# The Poincare series for Z/2Z summands in the homotopy groups of ko<2> (looped
# twice).
def poincare_series_pi_ko2_Z2():
  return one_m_tk_inv(8) + (t^7 * one_m_tk_inv(8))

# The Poincare series for Z summands in the homotopy groups of MSpin.
def poincare_series_pi_MSpin_Z(ko, ko2):
  return ((ko * poincare_series_pi_ko_Z())
    + (ko2 * poincare_series_pi_ko2_Z()))

# The Poincare series for Z/2Z summands in the homotopy groups of MSpin.
def poincare_series_pi_MSpin_Z2(ko, ko2):
  return ((ko * poincare_series_pi_ko_Z2())
    + (ko2 * poincare_series_pi_ko2_Z2())
    + poincare_series_MSpin_HZ2(ko, ko2))

# The Poincare series for Z summands in the homotopy groups of ku.
def poincare_series_pi_ku_Z():
  return one_m_tk_inv(2)

# The Poincare series for Z summands in the homotopy groups of MSpin^c.
def poincare_series_pi_MSpinc_Z(ku):
  return (ku * poincare_series_pi_ku_Z())

# The Poincare series for Z/2Z summands in the homotopy groups of MSpin^c.
def poincare_series_pi_MSpinc_Z2(ku):
  return poincare_series_MSpinc_HZ2(ku)

# The Poincare series for Z summands in the homotopy groups of ksp.
def poincare_series_pi_ksp_Z():
  return one_m_tk_inv(4)

# The Poincare series for Z/2Z summands in the homotopy groups of ksp.
def poincare_series_pi_ksp_Z2():
  return (t^5 * one_m_tk_inv(8)) + (t^6 * one_m_tk_inv(8))

# The homotopy groups of F are (abstractly) isomorphic to those of ko.

# The Poincare series for Z summands in the homotopy groups of MSpin^h.
def poincare_series_pi_MSpinh_Z(ksp, F):
  return ((ksp * poincare_series_pi_ksp_Z())
    + (F * poincare_series_pi_ko_Z()))

# The Poincare series for Z/2Z summands in the homotopy groups of MSpin^h.
def poincare_series_pi_MSpinh_Z2(ksp, F):
  return ((ksp * poincare_series_pi_ksp_Z2())
    + (F * poincare_series_pi_ko_Z2())
    + poincare_series_MSpinh_HZ2(ksp, F))


# --- Perform the calculation ---
if variant == "r":
  ko = poincare_series_MSpin_ko()
  ko2 = poincare_series_MSpin_ko2()
  Z = poincare_series_pi_MSpin_Z(ko, ko2)
  Z2 = poincare_series_pi_MSpin_Z2(ko, ko2)
  pi = []
  for i in range(0, n):
    pi.append([Z[i], Z2[i]])
  with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(pi)
elif variant == "c":
  ku = poincare_series_MSpinc_ku()
  Z = poincare_series_pi_MSpinc_Z(ku)
  Z2 = poincare_series_pi_MSpinc_Z2(ku)
  pi = []
  for i in range(0, n):
    pi.append([Z[i], Z2[i]])
  with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(pi)
elif variant == "h":
  ksp = poincare_series_MSpinh_ksp()
  F = poincare_series_MSpinh_F()
  Z = poincare_series_pi_MSpinh_Z(ksp, F)
  Z2 = poincare_series_pi_MSpinh_Z2(ksp, F)
  pi = []
  for i in range(0, n):
    pi.append([Z[i], Z2[i]])
  with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(pi)
