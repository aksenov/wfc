(ns wfc.core
  "Wave function collapse"
  (:require [clojure.set :as set]
            [clojure.spec.alpha :as s]
            [clojure.math :as math]))


(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)



;; sample loading
;; TODO in another module

;; sample analysis

(s/def :wfc/tile keyword?)
(s/def :wfc/sample-row (s/coll-of :wfc/tile :kind vector?))
(s/def :wfc/sample (s/coll-of :wfc/sample-row :kind vector?))


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


(defn entropy
  "Shannon's entropy"
  [variants weights]
  (let [tile-weights (map #(double (get weights %)) variants)
        tile-weight-sum (double (reduce + tile-weights))]
    (-
     (math/log tile-weight-sum)
     (/ (double (reduce + (map (fn [^double x] (* x (math/log x))) tile-weights)))
        tile-weight-sum))))

(defrecord Tile [r c collapsed? variants entropy index])

;; TODO get weights, variants from world and create new tile as world tile
(defn new-tile 
  "Create new tile"
  [r c variants weights]
  (->Tile r c false variants (entropy variants weights) [r c])
  #_{:r r :c c :collapsed? false :variants variants :entropy (entropy variants weights)})

(defn collapse-tile
  "Collapse tile to certain variant"
  [tile variant]
  (if (:collapsed? tile)
    tile
    (assoc tile
           :variants variant
           :collapsed? true
           :entropy nil)))


(defn get-tile
  "Get tile from the `world` in position [`r` `c`]"
  ([world [r c]]
   (get-tile world r c))
  ([world r c]
   (get-in world [:state r c])))


(defn lowest-entropy-tile
  "Search for the tile with lowest entropy"
  [world]
  (if (:lowest-entropy world)
    (assoc (get-tile world (:lowest-entropy-coords world))
           :index (:lowest-entropy-coords world))
    (reduce
     (fn [res el]
       (if (or (nil? res) (< ^double (:entropy el) ^double (:entropy res)))
         el
         res))
     nil
     (filter
      :entropy
      (flatten
       (map-indexed
        (fn [i row]
          (map-indexed (fn [j el]
                         (assoc el :index [i j])) row))
        (:state world)))))))

(defn lowest-entropy-tile
  "Search for the tile with lowest entropy"
  [world]
  (if (:lowest-entropy world)
    (get-tile world (:lowest-entropy-coords world))
    #_(assoc (get-tile world (:lowest-entropy-coords world))
           :index (:lowest-entropy-coords world))
    (reduce
     (fn [res el]
       (if (or (nil? res) (< ^double (:entropy el) ^double (:entropy res)))
         el
         res))
     nil
     (filter
      :entropy
      (flatten
       (:state world)
       #_(map-indexed
        (fn [i row]
          (map-indexed (fn [j el]
                         (assoc el :index [i j])) row))
        (:state world)))))))


(defn new-world
  [r c sample-analysis]
  (let [{:keys [adjacency weights superposition]} sample-analysis
        tile (new-tile 0 0 superposition weights)]
    {:state (into [] (for [r' (range r)]
                       (into [] (for [c' (range c)]
                                  (new-tile r' c' superposition weights))))
                  #_(repeat r (into [] (repeat c tile))))
     :rows r
     :cols c
     :collapsed? false
     :adjacency adjacency
     :weights weights
     :lowest-entropy (entropy superposition weights)
     :lowest-entropy-coords [(rand-int r) (rand-int c)]
     :superposition superposition}))


(defn adjacent-pos
  "Get coordinates of adjacents tiles left, right, up and down."
  [world ^long r ^long c]
  (let [left (dec c)
        right (inc c)
        up (dec r)
        down (inc r)
        {:keys [^long rows ^long cols]} world]
    {:left (when (>= left 0) [r left])
     :right (when (< right cols) [r right])
     :up (when (>= up 0) [up c])
     :down (when (< down rows) [down c])
     :left-up (when (and (>= up 0) (>= left 0)) [up left])
     :left-down (when (and (< down rows) (>= left 0)) [down left])
     :right-up (when (and (>= up 0) (< right cols)) [up right])
     :right-down (when (and (< down rows) (< right cols)) [down right])}))


(defn update-lowest-entropy
  "Remember lowest entropy and tile with it."
  ([world [r c] force?]
   (update-lowest-entropy world r c force?))
  ([world r c force?]
   (let [current-tile (get-tile world r c)
         ^double tile-entropy (:entropy current-tile)
         ^double lowest-entropy (:lowest-entropy world)]
     ;(prn "ULE" lowest-entropy tile-entropy (< tile-entropy lowest-entropy))
     (cond
       force? (assoc world
                     :lowest-entropy tile-entropy
                     :lowest-entropy-coords [r c])
       (nil? tile-entropy) world
       (and (nil? lowest-entropy) (some? tile-entropy))
       (assoc world
              :lowest-entropy tile-entropy
              :lowest-entropy-coords [r c])
       (< tile-entropy lowest-entropy) (assoc world
                                              :lowest-entropy tile-entropy
                                              :lowest-entropy-coords [r c])
       :else world)
     #_(if (or (<= lowest-entropy 0)
               (< tile-entropy lowest-entropy))
         (assoc world
                :lowest-entropy tile-entropy
                :lowest-entropy-coords [r c])
         world))))


(defn rand-weighted-choice
  "Random key from weight map."
  [weights]
  (let [wsum (reduce + (vals weights))
        r (long (rand-int wsum))]
    (loop [w (seq weights)
           s 0]
      (let [[k ^long v] (first w)]
        (if (or (empty? w)
                (<= r (+ v s)))
          k
          (recur (rest w)
                 (+ v s)))))))


(defn random-weighted-variant
  "Select variant from possible `variants` according to weight."
  [variants weights]
  ;(rand-nth (take 3 (reverse (sort-by #(get weights %) (seq (:variants tile))))))
  #_(rand-nth (seq (:variants tile)))
  (rand-weighted-choice (select-keys weights variants)))


(defn reduce-variants
  "Reduce variability of the tile by given variant"
  [tile variant weights]
  (let [{:keys [variants collapsed?]} tile]
    (if (or collapsed? (= 1 (count variants)))
      tile
      (let [new-variants (set/intersection variants variant)]
        (if (empty? new-variants)
          (do
            (println "-UB-" tile variants variant)
            (throw (ex-info "AAAA" {}))
            (assoc tile :uncertain? true))
          (assoc tile
                 :variants new-variants
                 :entropy (entropy new-variants weights)))))))

(defn adjacency-of-variants
  [variants adjacency]
  (reduce 
   (fn [res v]
     (merge-with set/union res (get adjacency v)))
   {:left #{}
    :right #{}
    :up #{}
    :down #{}}
   variants))


;; (defn collapse 
;;   "Collapse cell and reduce varability of it's adjacent tiles."
;;   [world r c variant]
;;   (let [{:keys [weights adjacency]} world
;;         adj-pos (adjacent-pos world r c)
;;         adj-tiles (get adjacency variant)
;;         slightly-collapsed-world (-> world
;;                                      (update-in [:state r c] collapse-tile variant)
;;                                      (update-lowest-entropy r c true))]
;;     (try
;;       (reduce
;;        (fn [w dir]
;;          (if-let [dir-coords (get adj-pos dir)]
;;            (let [[dr dc] dir-coords]
;;              (-> w
;;                  (update-in [:state dr dc] reduce-variants (get adj-tiles dir) weights)
;;                  (update-lowest-entropy dr dc false)))
;;            w))
;;        slightly-collapsed-world
;;        (shuffle [:left :right :up :down]))
;;       (catch Exception e))))

(defn update-tile-variants
  [world r c tile update-fn]
  (let [weights (:weights world)]
    (update-in world [:state r c] reduce-variants tile weights)))

(defn collapse
  "Collapse cell and reduce varability of it's adjacent tiles."
  [world r c variants]
  (loop [vs variants]
    (if (empty? vs)
      world ;; FIXME too permissive
      (let [variant (random-weighted-variant vs (:weights world))
            {:keys [weights adjacency]} world
            adj-pos (adjacent-pos world r c)
            adj-tiles (get adjacency variant)
            slightly-collapsed-world (-> world
                                         (update-in [:state r c] collapse-tile variant)
                                         (update-lowest-entropy r c true))
            world'
            (try
              (reduce
               (fn [w dir]
                 (if-let [dir-coords (get adj-pos dir)]
                   (let [[dr dc] dir-coords]
                     (-> w
                         (update-in [:state dr dc] reduce-variants (get adj-tiles dir) weights)
                         (update-lowest-entropy dr dc false)))
                   w))
               slightly-collapsed-world
               (shuffle [:left :right :up :down
                         ;:left-up :left-down :right-up :right-down
                         ]))
              (catch Exception e))]
        (if world'
          world'
          (recur (set/difference vs #{variant})))))))

(defn zero-entropy?
  [tile]
  (when tile
    (when (set? (:variants tile))
      (= 1 (count (:variants tile))))
    #_(when-let [^double e (:entropy tile)]
        (<= e 0))))

;; TODO progressive world weight

(defn collapse-adjacent
  ([world [r c] tile-variants]
   (collapse-adjacent world r c tile-variants))
  ([world r c tile-variants]
   (let [;v (if (set? variant) (first variant) variant)
         world'  (collapse world r c tile-variants)
         #_(loop [vs tile-variants]
                  ;(println "====> " [r c] vs)
                  (if (empty? vs)
                    ;(collapse world r c nil)
                    ;(throw (ex-info "WORLD FAILED" {}))
                    (do ;(println "=============" [r c] vs)
                        world)
                    (let [v (select-variant vs (:weights world))
                          www (collapse world r c v)]
                      (if www
                        www
                        (recur (set/difference vs #{v}))))))
         ;world' (collapse world r c v)
         {:keys [left right up down
                 left-up left-down right-up right-down]} (adjacent-pos world' r c)]
     (reduce
      (fn [w dir-coords]
        (let [dir-tile (when dir-coords (get-tile w dir-coords))]
          (if (zero-entropy? dir-tile)
            (do
              ;(prn "---" dir-coords (:variants dir-tile))
              (collapse-adjacent w dir-coords (:variants dir-tile)))
            (do
              ;(prn "+++" dir-coords (:variants dir-tile))
              w))))
      world'
      (shuffle [left right up down left-up left-down right-up right-down])))))


(defn collapse-world
  ([world]
   (collapse-world world (constantly true)))
  ([world continue-fn]
   (loop [world' world
          last-tile nil]
     (if (continue-fn world')
       (let [le-tile (lowest-entropy-tile world')]
         (if-not (or (nil? le-tile) (= le-tile last-tile))
           (recur (collapse-adjacent world' (:index le-tile) (:variants le-tile))
                  le-tile)
           (do (println le-tile)
               world')))
       world'))))


(defn print-collapsed
  [tile]
  (if (:uncertain? tile)
    (print (str (name :🟥)))
    (print (str (name (:variants tile))))))


(defn print-variants
  [tile]
  ;(print (str " " (count (:variants tile)) " "))
  (let [c (count (:variants tile))]
        (print (case c
                 0 "0️⃣"
                 1 "1️⃣"
                 2 "2️⃣"
                 3 "3️⃣"
                 4 "4️⃣"
                 5 "5️⃣"
                 6 "6️⃣"
                 "⬜"))))


(defn print-tile
  [tile]
  (if (:uncertain? tile)
    (print "❌")
    (if (:collapsed? tile)
      (print-collapsed tile)
      (print-variants tile))))


(defn print-world
  [world]
  (doseq [r (range (:rows world))]
    (doseq [c (range (:cols world))]
      (print-tile (get-tile world r c)))
    (println ":")))


(defn entropy
  "Shannon's entropy"
  [variants weights]
  ;(prn variants weights)
  (let [tile-weights (map #(get weights %) variants)
        ^double tile-weight-sum (reduce + tile-weights)]
    (-
     (Math/log tile-weight-sum)
     (/ (double (reduce + (map (fn [^double x] (* x (Math/log x))) tile-weights)))
        tile-weight-sum))))

(let [{:keys [weights]}
      (analyze-sample [[:🟫 :🟫 :🟫 :🟫]
                       [:🟩 :🟫 :🟫 :🟩]
                       [:🟦 :🟩 :🟩 :🟦]
                       [:🟦 :🟦 :🟦 :🟦]])]
  (entropy [:🟦] weights))
(comment
  (time
   (let [s1 [[:🟩 :🌳 :🟩 :🟩 :🌳  :🟡 :🟡 :🟩]
             [:🟩 :🟩 :🟩 :🟩 :🟡  :🟦 :🟨 :🟩]
             [:🌳 :🌳 :🌳 :🟡 :🟦  :🟦 :🟨 :🌳]
             [:🌲 :🌲 :🌳 :🟡 :🟦  :🟦 :🏨 :🌳]
             [:🌳 :🌲 :🟩 :🟩 :🟡  :🟦 :🟨 :🌲]
             [:🌳 :🌳 :🟩 :🟡 :🟦  :🟦 :🟨 :🌳]
             [:🌾 :🌳 :🟩 :🟡 :🟦  :🟨 :🌾 :🟩]
             [:🟩 :🌳 :🏡 :🟡 :🟨  :🟩 :🟩 :🟩]
             [:🟩 :🟩 :🌳 :🟩 :🟩  :🟩 :🟩 :🟩]]
         s2 [[:🟫 :🟫 :🟫 :🟫]
             [:🟩 :🟫 :🟫 :🟩]
             [:🟦 :🟩 :🟩 :🟦]
             [:🟦 :🟦 :🟦 :🟦]]
         s3 [[:🌾 :🌾]
             [:🌾 :🌾]]
         s4 [[:🟩 :🟩 :🟩 :🟩]
             [:🟩 :🟦 :🟦 :🟦]
             [:🟩 :🟦 :🌴 :🟦]
             [:🟩 :🟦 :🟦 :🟦]]
         s6 [[:🌳 :🌳 :🌳 :🌳 ]
             [:🌳 :🌲 :🌲 :🌳 ]
             [:🌳 :🌲 :🌲 :🌳 ]
             [:🌳 :🌳 :🌳 :🌳 ]]
          s7 [[:🟦 :🟦 :🟦]
              [:🟦 :🏝 :🟦]
              [:🟦 :🟦 :🟦]]
                   s8 [[:🟦 :🟦 :🟦]
                       [:🟦 :🟨 :🟦]
                       [:🟦 :🟦 :🟦]]
         s5 [[:A :A :A :A]
             [:A :A :A :A]
             [:A :A :A :A]
             [:A :C :C :A]
             [:C :B :B :C]
             [:C :B :B :C]
             [:C :B :B :C]]
         w (new-world 45 45 (analyze-samples [s1 s6 s7]))]
     (-> w
         ;(collapse-adjacent 0 0 [:🟫])
         ;(collapse-adjacent 1 1 [:🟨])
         ;(collapse-adjacent 2 2 [:🟦])
       ;(update-lowest-entropy 3 3)
         (collapse-world #(do
                            ;(println (:lowest-entropy-coords %))
                            ;(print-world %)
                            %

                            true))
      ;(collapse-adjacent 0 3 :🏡)
      ;(lowest-entropy-tile)
        (print-world)
         (println)
         )
  ;(collapse-adjacent w 1 1 :🟩)
  ;(lowest-entropy-tile )
     )))

;; TODO select with adjacency weights
;; TODO indexed search - X - use indexes inside tiles