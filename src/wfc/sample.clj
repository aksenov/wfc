(ns wfc.sample
  "Analyze sample."
  (:require [clojure.set :as set]))


(defn- calculate-adjacencies
  "Gather adjacency for tiles from `sample` for `direction1` and `direction2`."
  [sample direction1 direction2]
  (reduce
   (fn [acc row]
     (reduce
      (fn [acc2 [el er]]
        (-> acc2
            (update-in [el direction1] set/union #{er})
            (update-in [er direction2] set/union #{el})))
      acc
      (partition 2 1 row)))
   {}
   sample))


(defn- transpose
  "Transpose `sample` matrix."
  [sample]
  (apply (partial mapv vector) sample))


(defn adjacency 
  "Make a map of tile to left, right, up and down adjacent tile sets. "
  [sample]
  (merge-with
   (partial merge-with conj)
   (calculate-adjacencies sample :right :left)
   (calculate-adjacencies (transpose sample) :down :up)))


(defn weights
  "Specific weight of the tile in the sample"
  [sample]
  (frequencies (flatten sample)))


(defn superposition
  "All possible sample tiles."
  [sample]
  (into #{} (flatten sample)))


(defn analyze-sample
  "Calculate sample's main parameters: adjacency, weight and superposition."
  [sample]
  {:adjacency (adjacency sample)
   :weights (weights sample)
   :superposition (superposition sample)})



(defn analyze-samples
  "Do analyze multiple samples and joins results."
  [samples]
  (reduce
   (fn [res s]
     (-> res
         (update :adjacency #(merge-with
                              (fn [r e] (merge-with into r e)) % (:adjacency s)))
         (update :superposition #(into % (:superposition s)))
         (update :weights #(merge-with + % (:weights s)))))
   {:superposition #{}
    :adjacency {}
    :weights {}}
   (map analyze-sample samples)))
